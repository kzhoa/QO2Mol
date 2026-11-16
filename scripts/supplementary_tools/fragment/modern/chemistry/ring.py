# Shared structural helper extracted from the original implementation.
from collections import OrderedDict
from copy import deepcopy
from .utils import combine_arr
from .aromatics import Aromatics


def find_next_nodel(a, connect_dict: OrderedDict, nodel, path_arr):
    connect_atom = list(set(connect_dict[a]) - set(path_arr[1:]))
    for b in connect_atom:
        path_arr.append(b)
        if b not in nodel:
            path_arr = find_next_nodel(b, connect_dict, nodel, path_arr)
    return path_arr


def find_cyclo(connect_dict: OrderedDict) -> list:
    if not connect_dict:
        return []
    # Cycle searches operate on connected graphs. A union of independent
    # degree-two components must not become a single fictitious ring.
    unseen = set(connect_dict)
    components = []
    for start in sorted(connect_dict):
        if start not in unseen:
            continue
        unseen.remove(start)
        pending, component = [start], []
        while pending:
            a = pending.pop()
            component.append(a)
            for b in connect_dict[a]:
                if b in unseen:
                    unseen.remove(b)
                    pending.append(b)
        components.append(component)
    if len(components) > 1:
        return [
            cycle
            for component in components
            for cycle in find_cyclo({a: connect_dict[a] for a in sorted(component)})
        ]
    if all(len(neighbors) == 2 for neighbors in connect_dict.values()):
        start = min(connect_dict)
        path = [start]
        previous, current = start, min(connect_dict[start])
        while current != start:
            path.append(current)
            following = next(a for a in connect_dict[current] if a != previous)
            previous, current = current, following
        return [path + [start]]
    double, nodel, edge, k, ring_arr, tmp_ring_arr = [], [], 0, 0, [], []
    edge = sum([len(ele) for ele in connect_dict.values()])
    double = [key for key in connect_dict.keys() if len(connect_dict[key]) == 2]  # noqa
    nodel = [key for key in connect_dict.keys() if len(connect_dict[key]) != 2]
    face = edge // 2 + 2 - len(connect_dict)
    for a in nodel:
        this_face = len(connect_dict[a])
        m, path_tmp, same_tmp, pp = 0, {}, [], 0
        for i in range(this_face):
            aaa = connect_dict[a][i]
            for j in range(len(connect_dict[aaa])):
                if connect_dict[aaa][j] != a:
                    path_tmp[str(pp)] = [a, aaa, connect_dict[aaa][j]]
                    pp += 1
        for i in range(20):
            path_tmp1 = deepcopy(path_tmp)
            for n, patha in path_tmp1.items():
                connect_atom = list(set(connect_dict[patha[-1]]) - set(patha[1:]))
                for j in range(len(connect_atom)):
                    b = connect_atom[j]
                    if j != 0:
                        pp += 1
                    x = n if j == 0 else str(pp)
                    tmp_arr = deepcopy(patha) + [b]
                    if b not in nodel:
                        tmp_arr = find_next_nodel(b, connect_dict, nodel, tmp_arr)
                    if tmp_arr[-1] == a:
                        equal_set, is_subset = False, False
                        for aaa in same_tmp:
                            if set(tmp_arr) == set(aaa):
                                equal_set = True
                                break
                            elif set(tmp_arr) > set(aaa):
                                is_subset = True
                                break
                        if not (equal_set or is_subset):
                            same_tmp.append(tmp_arr)
                            m += 1
                            equal_set = False
                            for aaa in ring_arr:
                                if set(aaa) == set(tmp_arr):
                                    equal_set = True
                                    break
                            if not equal_set:
                                k += 1
                                ring_arr.append(tmp_arr)
                            if x == n:
                                path_tmp.pop(x)
                        elif not equal_set and is_subset:
                            m += 1
                    else:
                        path_tmp[x] = tmp_arr
            while True:
                this_flag, break_flag = False, True
                for i in range(len(ring_arr)):
                    for j in range(len(ring_arr)):
                        if i == j:
                            continue
                        if set(nodel) - set(ring_arr[i]) == set(nodel) - set(
                            ring_arr[j]
                        ):
                            if len(ring_arr[i]) > len(ring_arr[j]):
                                tmp_ring_arr.append(ring_arr[i])
                                del ring_arr[i]
                            else:
                                tmp_ring_arr.append(ring_arr[j])
                                del ring_arr[j]
                            this_flag = True
                            break_flag = False
                            break
                    if this_flag:
                        break
                if break_flag:
                    break
            if m >= this_face:
                break
    if face - len(ring_arr) >= 2:
        aa1 = set([ele for sublist in ring_arr for ele in sublist])
        aa2 = set([ele for sublist in ring_arr + tmp_ring_arr for ele in sublist])
        diff_set = aa2 - aa1
        lose_total = []
        for a in diff_set:
            if a in lose_total:
                continue
            for rr in tmp_ring_arr:
                if a in rr:
                    ring_arr.append(rr)
                    lose_total.extend(rr)
                    break
    return ring_arr


def cyclo_property(molecule, cyclos):
    total_ring_arr = []
    for arr in cyclos:
        total_ring_arr.extend(arr)
    for arr in cyclos:
        aromobj = Aromatics(molecule, arr, total_ring_arr)
        arom = aromobj.get_aromatics()
        arr[-1] = arom
    return cyclos


def cyclo_blocks(cyclos, connectivity=None, ring_atoms=None):
    rings = {}
    # ring_stru = []
    for arr in cyclos:
        n = len(arr) - 1
        nn = int(n / 2)
        ss = f"R{n}-{arr[0]}_{arr[nn]}_{arr[-1]}"
        suffix = 0
        while ss in rings:
            suffix += 1
            ss = f"R{n}-{arr[0]}_{arr[nn]}--{suffix}_{arr[-1]}"
        rings[ss] = arr

    if connectivity is not None:
        used_rings = deepcopy(rings)
        for aa, bb in rings.items():
            tmp_arr = deepcopy(bb[:-1])
            for ii in bb[:-1]:
                tmp_arr.extend(
                    [jj for jj in connectivity[ii] if jj in ring_atoms and jj not in bb]
                )
            tmp_arr.append(bb[-1])
            used_rings[aa] = tmp_arr
    else:
        used_rings = rings

    tmp = []
    for aa, bb in used_rings.items():
        tmp.append(bb[:-1])
    ring_blocks = combine_arr(tmp)
    ring_block_components = [[] for i in range(len(ring_blocks))]
    for aa, bb in rings.items():
        for ii, rr in enumerate(ring_blocks):
            if set(bb[:-1]).issubset(set(rr)):
                ring_block_components[ii].append(aa)
                break
    return rings, ring_blocks, ring_block_components
