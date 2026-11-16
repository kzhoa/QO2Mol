"""Non-periodic molecular graphs. Edges point from source j to target i."""
from importlib.util import find_spec

import torch
from qqtools.torch.qmgraph import qtriplets
from qqtools.torch.qscatter import scatter as qq_scatter


@torch.no_grad()
def glorot_orthogonal(tensor, scale):
    torch.nn.init.orthogonal_(tensor)
    tensor.mul_((scale / ((tensor.size(-2) + tensor.size(-1)) * tensor.var())).sqrt())


def scatter(values, index, dim=0, dim_size=None, reduce="sum"):
    if index.numel() == 0:
        shape = list(values.shape)
        shape[dim] = dim_size or 0
        return values.sum(dim, keepdim=True).expand(shape)
    return qq_scatter(values, index, dim, dim_size, reduce)


def expand_ranges(starts, lengths):
    """Return (group index, item index) for consecutive ranges."""
    groups = torch.arange(len(starts), device=starts.device).repeat_interleave(lengths)
    offsets = lengths.cumsum(0) - lengths
    local = torch.arange(groups.numel(), device=starts.device) - offsets[groups]
    return groups, starts[groups] + local


@torch.no_grad()
def radius_graph(pos, r, batch=None, max_num_neighbors=32):
    """Keep the first in-radius sources by node index, excluding the target."""
    cap = min(max_num_neighbors, len(pos) - 1)
    if cap <= 0:
        return torch.empty((2, 0), device=pos.device, dtype=torch.long)
    if batch is None:
        batch = torch.zeros(len(pos), dtype=torch.long, device=pos.device)
    # Large caps can make a padded neighbor buffer larger than the pair list.
    if (
        pos.is_cuda
        and cap <= 256
        and pos.shape[1] == 3
        and pos.dtype in (torch.float32, torch.float64)
        and torch.cuda.get_device_capability(pos.device)[0] >= 7
        and find_spec("triton") is not None
    ):
        from .cuda_neighbors import cuda_radius_graph

        return cuda_radius_graph(pos, r, batch.contiguous(), cap)
    order = torch.argsort(batch, stable=True)
    counts = batch.bincount()
    starts = counts.cumsum(0) - counts
    target, source = expand_ranges(starts[batch], counts[batch])
    source = order[source]
    keep = (pos[source] - pos[target]).square().sum(-1) < r * r
    keep &= source != target
    source, target = source[keep], target[keep]
    degree = target.bincount(minlength=len(pos))
    first = (degree.cumsum(0) - degree)[target]
    keep = torch.arange(len(target), device=target.device) - first < cap
    return torch.stack((source[keep], target[keep]))


def triplets(edge_index, num_nodes):
    offsets = edge_index.new_zeros((edge_index.size(1), 3))
    return qtriplets(edge_index, offsets, num_nodes)
