"""Bounded-memory CUDA radius search for three-dimensional molecules."""
import torch
import triton
import triton.language as tl


@triton.jit
def _radius_kernel(pos, order, starts, stops, neighbors, radius_sq: tl.constexpr,
                   cap: tl.constexpr, block: tl.constexpr):
    target = tl.program_id(0)
    start = tl.load(starts + target)
    stop = tl.load(stops + target)
    x = tl.load(pos + target * 3)
    y = tl.load(pos + target * 3 + 1)
    z = tl.load(pos + target * 3 + 2)
    cutoff = tl.full((), radius_sq, x.dtype)
    lanes = tl.arange(0, block)
    count = 0
    while (start < stop) & (count < cap):
        valid = start + lanes < stop
        source = tl.load(order + start + lanes, mask=valid, other=0)
        dx = tl.load(pos + source * 3) - x
        dy = tl.load(pos + source * 3 + 1) - y
        dz = tl.load(pos + source * 3 + 2) - z
        inside = valid & (source != target) & (dx * dx + dy * dy + dz * dz < cutoff)
        rank = count + tl.cumsum(inside.to(tl.int32), 0) - 1
        tl.store(neighbors + target * cap + rank, source, mask=inside & (rank < cap))
        count += tl.sum(inside.to(tl.int32), 0)
        start += block


def cuda_radius_graph(pos, r, batch, cap):
    graph, order = torch.sort(batch, stable=True)
    starts = torch.searchsorted(graph, batch)
    stops = torch.searchsorted(graph, batch, right=True)
    neighbors = order.new_full((len(pos), cap), -1)
    with torch.cuda.device(pos.device):
        _radius_kernel[(len(pos),)](pos.contiguous(), order, starts, stops,
                                    neighbors, r * r, cap, 64, num_warps=2,
                                    enable_fp_fusion=False)
    target, slot = (neighbors >= 0).nonzero(as_tuple=True)
    return torch.stack((neighbors[target, slot], target))
