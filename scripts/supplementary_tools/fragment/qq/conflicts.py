"""Resolve local merge conflicts only when the admissible partition is unique."""

from copy import deepcopy
from dataclasses import replace
from itertools import combinations
from rdkit import Chem
from .basic_units import partition_basic_units as baseline
from .partition import _baseline_partition
from .models import ElementaryFragment

LIMIT = 16


def partition_key(groups):
    return tuple(sorted((tuple(sorted(g)) for g in groups)))


def merged_groups(blocks, unions):
    owner = {i: i for i in range(len(blocks))}
    for support in unions:
        labels = {owner[i] for i in support}
        target = min(labels)
        for i in owner:
            if owner[i] in labels:
                owner[i] = target
    result = {}
    for i, atoms in enumerate(blocks):
        result.setdefault(owner[i], set()).update(atoms)
    return partition_key(result.values())


def metrics(groups, heavy, edges):
    owner = {a: i for i, g in enumerate(groups) for a in g}
    ports = [0] * len(groups)
    for a, b in edges:
        if owner[a] != owner[b]:
            ports[owner[a]] += 1
            ports[owner[b]] += 1
    sizes = [sum((a in heavy for a in g)) for g in groups]
    return (
        sum((0 < n <= 2 and p >= 2 for n, p in zip(sizes, ports))),
        sum((0 < n <= 2 for n in sizes)),
        sum((16 <= n <= 20 for n in sizes)),
        max(sizes),
    )


def local_metrics(groups, heavy, edges, scope):
    owner = {a: i for i, g in enumerate(groups) for a in g}
    ports = [0] * len(groups)
    for a, b in edges:
        if owner[a] != owner[b]:
            ports[owner[a]] += 1
            ports[owner[b]] += 1
    chosen = [
        (sum((a in heavy for a in g)), ports[i])
        for i, g in enumerate(groups)
        if set(g) & scope
    ]
    assert all((set(g) <= scope for g in groups if set(g) & scope))
    return (
        sum((0 < n <= 2 and p >= 2 for n, p in chosen)),
        sum((0 < n <= 2 for n, p in chosen)),
        sum((16 <= n <= 20 for n, p in chosen)),
        max((n for n, p in chosen)),
    )


def components(proposals):
    remaining = {i for i, p in enumerate(proposals) if p["conflicts"]}
    result = []
    while remaining:
        group = {min(remaining)}
        todo = list(group)
        while todo:
            for j in proposals[todo.pop()]["conflicts"]:
                if j not in group:
                    group.add(j)
                    todo.append(j)
        remaining -= group
        result.append(sorted(group))
    return result


def solve_group(blocks, heavy, edges, proposals, indices, fixed, reference):
    scope = set().union(
        *(blocks[j] for i in indices for j in proposals[i]["baseline_blocks"])
    )
    reference = local_metrics(merged_groups(blocks, fixed), heavy, edges, scope)
    log = dict(
        scope=sorted(scope),
        reference=list(reference),
        indices=indices,
        node_count=len(indices),
        eligible=[i for i in indices if proposals[i]["heavy_atoms"] <= 20],
        schemes=[],
        accepted_unions=[],
    )
    if len(indices) > LIMIT:
        log.update(
            status="defer_search_limit", compatible_schemes=0, frontier_partitions=0
        )
        return log
    eligible = log["eligible"]
    accepted = {}
    for count in range(len(eligible) + 1):
        for subset in combinations(eligible, count):
            used = set()
            valid = True
            for i in subset:
                support = set(proposals[i]["baseline_blocks"])
                if used & support:
                    valid = False
                    break
                used |= support
            if not valid:
                continue
            unions = [proposals[i]["baseline_blocks"] for i in subset]
            groups = merged_groups(blocks, [*fixed, *unions])
            m = local_metrics(groups, heavy, edges, scope)
            reasons = []
            if m[0] >= reference[0]:
                reasons.append("no_multiport_tiny_reduction")
            if m[1] > reference[1]:
                reasons.append("tiny_increase")
            if m[2] > reference[2]:
                reasons.append("extended_increase")
            if m[3] > reference[3]:
                reasons.append("local_max_increase")
            entry = dict(
                subset=list(subset),
                metrics=list(m),
                rejections=reasons,
                partition=[list(g) for g in groups],
            )
            log["schemes"].append(entry)
            if reasons:
                continue
            item = accepted.setdefault(groups, dict(metrics=m, subsets=[]))
            item["subsets"].append(list(subset))
    frontier = {
        g: x
        for g, x in accepted.items()
        if not any(
            (
                all((a <= b for a, b in zip(y["metrics"], x["metrics"])))
                and any((a < b for a, b in zip(y["metrics"], x["metrics"])))
                for y in accepted.values()
            )
        )
    }
    log["compatible_schemes"] = len(log["schemes"])
    log["frontier_partitions"] = len(frontier)
    log["frontier"] = [
        dict(
            partition=[list(a) for a in g],
            metrics=list(x["metrics"]),
            subsets=x["subsets"],
        )
        for g, x in frontier.items()
    ]
    if len(frontier) != 1:
        log["status"] = (
            "defer_no_admissible_scheme" if not frontier else "defer_ambiguous_frontier"
        )
        return log
    groups, item = next(iter(frontier.items()))
    original = {frozenset(g): i for i, g in enumerate(blocks)}
    fixed_groups = set(merged_groups(blocks, fixed))
    unions = []
    for group in groups:
        if group in fixed_groups:
            continue
        members = set(group)
        support = [i for i, g in enumerate(blocks) if set(g) <= members]
        assert (
            len(support) > 1
            and set().union(*(set(blocks[i]) for i in support)) == members
        )
        assert sum((a in heavy for a in members)) <= 20
        unions.append(support)
    log.update(
        status="applied_unique_partition",
        accepted_unions=unions,
        contributing_proposals=sorted(
            {i for subset in item["subsets"] for i in subset}
        ),
    )
    return log


def resolve(mol, config):
    old = baseline(mol, config)
    original, _, _, _ = _baseline_partition(mol, config)
    proposals = deepcopy(old[4])
    blocks = [set(g.atoms) for g in original]
    heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    fixed = [p["baseline_blocks"] for p in proposals if p["status"] == "applied"]
    reference = metrics([g.atoms for g in old[0]], heavy, edges)
    logs = [
        solve_group(blocks, heavy, edges, proposals, group, fixed, reference)
        for group in components(proposals)
    ]
    new_unions = [u for log in logs for u in log["accepted_unions"]]
    supports = [set(u) for u in new_unions]
    assert all((not a & b for i, a in enumerate(supports) for b in supports[i + 1 :]))
    target = merged_groups(blocks, [*fixed, *new_unions])
    final_metric = metrics(target, heavy, edges)
    assert all((a <= b for a, b in zip(final_metric, reference)))
    if not new_unions:
        return (old, logs)
    owner = {a: i for i, g in enumerate(target) for a in g}
    decisions = []
    new_cuts_cancelled = []
    for d in old[2]:
        if d.cut and owner[d.atoms[0]] == owner[d.atoms[1]]:
            new_cuts_cancelled.append(list(d.atoms))
            d = replace(
                d,
                cut=False,
                reason="basic_unit_conflict_resolved",
                protections=(*d.protections, "basic_unit_conflict_resolved"),
            )
        decisions.append(d)
    graph = Chem.RWMol(mol)
    for d in decisions:
        if d.cut:
            graph.RemoveBond(*d.atoms)
    groups = Chem.GetMolFrags(graph, asMols=False, sanitizeFrags=False)
    assert partition_key(groups) == target
    fragments = [
        ElementaryFragment(
            i,
            tuple(g),
            tuple((a for a in g if a in heavy)),
            any((mol.GetAtomWithIdx(a).IsInRing() for a in g)),
        )
        for i, g in enumerate(groups)
    ]
    owners = {a: g.index for g in fragments for a in g.atoms}
    adj = {g.index: set() for g in fragments}
    for d in decisions:
        if d.cut:
            a, b = [owners[x] for x in d.atoms]
            assert a != b
            adj[a].add(b)
            adj[b].add(a)
    for log in logs:
        if log["status"] == "applied_unique_partition":
            for i in log["contributing_proposals"]:
                p = proposals[i]
                p["status"] = "applied_conflict_resolved"
                p["cancelled_bonds"] = [
                    e for e in new_cuts_cancelled if set(e) <= set(p["members"])
                ]
    return ((fragments, adj, decisions, old[3], proposals), logs)
