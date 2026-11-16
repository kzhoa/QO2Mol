"""Bounded carbon periphery at aromatic ethers."""

from dataclasses import replace
from rdkit import Chem
from .models import ElementaryFragment

QUERY = Chem.MolFromSmarts("[c]-[O;+0;R0]-[C;+0;X4;R0]")


def endpoints(mol, core):
    out = {}
    for ar, o, c in mol.GetSubstructMatches(QUERY, maxMatches=100001):
        oxygen = mol.GetAtomWithIdx(o)
        if (
            len([n for n in oxygen.GetNeighbors() if n.GetAtomicNum() > 1]) != 2
            or oxygen.GetNumRadicalElectrons()
            or o in core
        ):
            continue
        out[o] = ar
    return out


def resolve(mol, config, previous):
    groups, _, decisions, motifs, operations = previous
    owner = {a: g.index for g in groups for a in g.atoms}
    byindex = {g.index: g for g in groups}
    core = {a for m in motifs for a in m["atoms"]}
    coreblocks = {owner[a] for a in core}
    ep = endpoints(mol, core)
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
        anchors = [(a, b) for a, b in boundary if b in ep]
        if not anchors:
            continue
        reasons = []
        if len(boundary) < 2:
            reasons.append("not_multiport")
        if members & core:
            reasons.append("contains_core")
        if any(
            (
                mol.GetAtomWithIdx(a).GetAtomicNum() != 6
                or mol.GetAtomWithIdx(a).GetFormalCharge()
                or mol.GetAtomWithIdx(a).GetNumRadicalElectrons()
                or mol.GetAtomWithIdx(a).IsInRing()
                or (
                    mol.GetAtomWithIdx(a).GetHybridization()
                    != Chem.HybridizationType.SP3
                )
                for a in g.heavy_atoms
            )
        ):
            reasons.append("not_neutral_acyclic_sp3_carbon")
        if len(anchors) != 1:
            reasons.append("multiple_aryl_ether_anchors")
        target = owner[anchors[0][1]] if len(anchors) == 1 else None
        if target is not None:
            o = anchors[0][1]
            if owner[ep[o]] != target:
                reasons.append("aryl_O_not_preserved")
            if sum((owner[b] == target for a, b in boundary)) != 1:
                reasons.append("multiple_target_edges")
            for a, b in boundary:
                if b == o:
                    continue
                atom = mol.GetAtomWithIdx(b)
                if (
                    atom.GetAtomicNum() != 6
                    or atom.GetFormalCharge()
                    or atom.GetNumRadicalElectrons()
                    or (owner[b] in coreblocks)
                ):
                    reasons.append("competing_hetero_or_core")
                    break
            if (members | set(byindex[target].atoms)) & unresolved:
                reasons.append("unresolved_merge_claim")
        item = dict(
            block_index=g.index,
            members=list(g.atoms),
            boundary=[list(e) for e in boundary],
            anchors=[list(e) for e in anchors],
            target=target,
            rejections=reasons,
        )
        items.append(item)
        if not reasons:
            eligible.setdefault(target, []).append(g.index)
    batches = []
    for target, donors in sorted(eligible.items()):
        support = [target, *sorted(donors)]
        members = set().union(*(set(byindex[i].atoms) for i in support))
        added = sum((len(byindex[i].heavy_atoms) for i in donors))
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
            )
        )
    final_owner = dict(owner)
    seen = set()
    for batch in batches:
        if batch["status"] != "applied":
            continue
        assert not seen & set(batch["members"])
        seen.update(batch["members"])
        for a in batch["members"]:
            final_owner[a] = batch["target"]
    ds = [
        (
            replace(
                d,
                cut=False,
                reason="aryl_ether_unique_periphery",
                protections=(*d.protections, "aryl_ether_unique_periphery"),
            )
            if d.cut and final_owner[d.atoms[0]] == final_owner[d.atoms[1]]
            else d
        )
        for d in decisions
    ]
    graph = Chem.RWMol(mol)
    for d in ds:
        if d.cut:
            graph.RemoveBond(*d.atoms)
    gs = [
        ElementaryFragment(
            i,
            tuple(g),
            tuple((a for a in g if mol.GetAtomWithIdx(a).GetAtomicNum() > 1)),
            any((mol.GetAtomWithIdx(a).IsInRing() for a in g)),
        )
        for i, g in enumerate(
            Chem.GetMolFrags(graph, asMols=False, sanitizeFrags=False)
        )
    ]
    final = {a: g.index for g in gs for a in g.atoms}
    adj = {g.index: set() for g in gs}
    for d in ds:
        if d.cut:
            a, b = [final[i] for i in d.atoms]
            assert a != b
            adj[a].add(b)
            adj[b].add(a)
    return (
        (gs, adj, ds, motifs, operations),
        dict(endpoints=ep, periphery=items, batches=batches),
    )
