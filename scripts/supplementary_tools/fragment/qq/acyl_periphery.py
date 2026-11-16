"""Bounded carbon periphery at acyl nitrogen and oxygen."""

from dataclasses import replace
from rdkit import Chem
from .models import ElementaryFragment

QUERY = Chem.MolFromSmarts("[C;+0](=[O;+0])-[N,O;+0]")


def templates(mol):
    matches = mol.GetSubstructMatches(QUERY, maxMatches=100001)
    by_center = {}
    for c, o, x in matches:
        by_center.setdefault(c, dict(oxygen=o, hetero=set()))["hetero"].add(x)
    result = []
    for c, data in by_center.items():
        hetero = sorted(data["hetero"])
        elements = sorted((mol.GetAtomWithIdx(i).GetAtomicNum() for i in hetero))
        kind = {(7,): "amide", (8,): "ester", (7, 7): "urea", (7, 8): "carbamate"}.get(
            tuple(elements)
        )
        if not kind or any((mol.GetAtomWithIdx(i).IsInRing() for i in [c, *hetero])):
            continue
        result.append(dict(kind=kind, carbonyl=c, oxygen=data["oxygen"], hetero=hetero))
    return result


def resolve(mol, config, previous):
    groups, _, decisions, motifs, operations = previous
    byindex = {g.index: g for g in groups}
    owner = {a: g.index for g in groups for a in g.atoms}
    all_core = {a for m in motifs for a in m["atoms"]}
    core_blocks = {owner[a] for a in all_core}
    choices = templates(mol)
    endpoints = {}
    for t in choices:
        for x in t["hetero"]:
            endpoints.setdefault(x, []).append(t)
    unresolved = {
        a for p in operations if p["status"].startswith("defer") for a in p["members"]
    }
    cuts = [d.atoms for d in decisions if d.cut]
    items = []
    eligible = {}
    for g in groups:
        if not 1 <= len(g.heavy_atoms) <= 2:
            continue
        members = set(g.atoms)
        boundary = [
            (a, b) if a in members else (b, a)
            for a, b in cuts
            if (a in members) != (b in members)
        ]
        neighbors = {owner[b] for a, b in boundary}
        near = neighbors & core_blocks
        if not near:
            continue
        reasons = []
        if len(boundary) < 2:
            reasons.append("not_multiport")
        if members & all_core:
            reasons.append("contains_core")
        if any(
            (
                mol.GetAtomWithIdx(a).GetAtomicNum() != 6
                or mol.GetAtomWithIdx(a).GetFormalCharge() != 0
                or mol.GetAtomWithIdx(a).IsInRing()
                or (
                    mol.GetAtomWithIdx(a).GetHybridization()
                    != Chem.HybridizationType.SP3
                )
                for a in g.heavy_atoms
            )
        ):
            reasons.append("not_neutral_acyclic_sp3_carbon")
        if len(near) != 1:
            reasons.append("multiple_core_neighbors")
        target = next(iter(near)) if len(near) == 1 else None
        edges = [(a, b) for a, b in boundary if owner[b] == target]
        if len(edges) != 1:
            reasons.append("not_single_anchor_edge")
        anchor = edges[0][1] if len(edges) == 1 else None
        if anchor not in endpoints:
            reasons.append("not_supported_acyl_NO_endpoint")
        if target is not None and (members | set(byindex[target].atoms)) & unresolved:
            reasons.append("unresolved_merge_claim")
        item = dict(
            block_index=g.index,
            members=list(g.atoms),
            heavy_atoms=len(g.heavy_atoms),
            boundary=[list(e) for e in boundary],
            core_neighbors=sorted(near),
            target=target,
            anchor=anchor,
            rejections=reasons,
        )
        if anchor in endpoints:
            item["templates"] = endpoints[anchor]
        items.append(item)
        if not reasons:
            eligible.setdefault(target, []).append(item)
    batches = []
    for target, linkers in sorted(eligible.items()):
        support = [target, *sorted((x["block_index"] for x in linkers))]
        members = set().union(*(set(byindex[i].atoms) for i in support))
        added = sum((x["heavy_atoms"] for x in linkers))
        heavy = sum((mol.GetAtomWithIdx(a).GetAtomicNum() > 1 for a in members))
        reasons = []
        if added > 4:
            reasons.append("extension_over4")
        if heavy > 15:
            reasons.append("union_over15")
        batches.append(
            dict(
                target=target,
                blocks=support,
                members=sorted(members),
                added_heavy=added,
                heavy_atoms=heavy,
                status="deferred" if reasons else "applied",
                rejections=reasons,
                linkers=[x["block_index"] for x in linkers],
            )
        )
    applied = [b for b in batches if b["status"] == "applied"]
    seen = set()
    final_owner = dict(owner)
    for batch in applied:
        assert not seen & set(batch["members"])
        seen.update(batch["members"])
        for a in batch["members"]:
            final_owner[a] = batch["target"]
    final_decisions = []
    for d in decisions:
        if d.cut and final_owner[d.atoms[0]] == final_owner[d.atoms[1]]:
            d = replace(
                d,
                cut=False,
                reason="acyl_unique_NO_periphery",
                protections=(*d.protections, "acyl_unique_NO_periphery"),
            )
        final_decisions.append(d)
    graph = Chem.RWMol(mol)
    for d in final_decisions:
        if d.cut:
            graph.RemoveBond(*d.atoms)
    components = Chem.GetMolFrags(graph, asMols=False, sanitizeFrags=False)
    expected = {}
    for a, label in final_owner.items():
        expected.setdefault(label, set()).add(a)
    assert {frozenset(g) for g in components} == {
        frozenset(g) for g in expected.values()
    }
    fragments = [
        ElementaryFragment(
            i,
            tuple(g),
            tuple((a for a in g if mol.GetAtomWithIdx(a).GetAtomicNum() > 1)),
            any((mol.GetAtomWithIdx(a).IsInRing() for a in g)),
        )
        for i, g in enumerate(components)
    ]
    final = {a: g.index for g in fragments for a in g.atoms}
    adj = {g.index: set() for g in fragments}
    for d in final_decisions:
        if d.cut:
            a, b = [final[i] for i in d.atoms]
            assert a != b
            adj[a].add(b)
            adj[b].add(a)
    return (
        (fragments, adj, final_decisions, motifs, operations),
        dict(templates=choices, periphery=items, batches=batches),
    )
