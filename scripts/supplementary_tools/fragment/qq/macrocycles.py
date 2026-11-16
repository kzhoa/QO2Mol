"""Unique paired cuts of isolated nonaromatic macrocycles."""

from dataclasses import replace
from itertools import combinations
from rdkit import Chem
from .models import ElementaryFragment
from .products import cap


def ring_systems(mol):
    adj = {a.GetIdx(): set() for a in mol.GetAtoms() if a.IsInRing()}
    for b in mol.GetBonds():
        if b.IsInRing():
            a, c = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            adj[a].add(c)
            adj[c].add(a)
    remaining = set(adj)
    out = []
    while remaining:
        todo = [min(remaining)]
        seen = set()
        while todo:
            a = todo.pop()
            if a in seen:
                continue
            seen.add(a)
            todo.extend(adj[a] - seen)
        remaining -= seen
        out.append(
            dict(
                atoms=sorted(seen),
                size=len(seen),
                simple_cycle=all((len(adj[a]) == 2 for a in seen)),
                aromatic_atoms=sum(
                    (mol.GetAtomWithIdx(a).GetIsAromatic() for a in seen)
                ),
                cycle_rank=sum((len(adj[a]) for a in seen)) // 2 - len(seen) + 1,
            )
        )
    return out


def partition_key(parts):
    return tuple(sorted((tuple(sorted(p)) for p in parts)))


def split(mol, members, existing, new):
    rw = Chem.RWMol(mol)
    for e in [*existing, *new]:
        rw.RemoveBond(*e)
    return [
        list(p)
        for p in Chem.GetMolFrags(rw, asMols=False, sanitizeFrags=False)
        if set(p) & set(members)
    ]


def stereo_present(product, construction):
    return {
        a
        for i, a in enumerate(construction["product_to_parent"])
        if a is not None
        and product.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    }


def resolve(mol, config, previous):
    groups, _, ds, motifs, operations = previous
    existing = [list(d.atoms) for d in ds if d.cut]
    core = {a for t in motifs for a in t["atoms"]}
    systems = ring_systems(mol)
    decisions = {frozenset(d.atoms): d for d in ds}
    merged = {
        a
        for d in ds
        if any(
            (
                p.startswith("basic_unit_")
                or p.endswith("_periphery")
                or p.startswith("conflict_")
                for p in d.protections
            )
        )
        for a in d.atoms
    }
    records = []
    chosen = []
    for g in groups:
        if len(g.heavy_atoms) <= 20:
            continue
        members = set(g.atoms)
        rings = [
            r
            for r in systems
            if set(r["atoms"]) <= members
            and r["simple_cycle"]
            and (r["size"] >= 12)
            and (not r["aromatic_atoms"])
        ]
        if not rings:
            continue
        sets = []
        all_edges = []
        for ring in rings:
            nodes = set(ring["atoms"])
            edges = []
            for bond in mol.GetBonds():
                a, b = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                d = decisions[frozenset((a, b))]
                if (
                    a not in nodes
                    or b not in nodes
                    or (not bond.IsInRing())
                    or d.cut
                    or (set(d.protections) != {"ring_bond"})
                    or (bond.GetBondType() != Chem.BondType.SINGLE)
                    or bond.GetIsConjugated()
                ):
                    continue
                if all(
                    (
                        mol.GetAtomWithIdx(x).GetAtomicNum() == 6
                        and (not mol.GetAtomWithIdx(x).GetFormalCharge())
                        and (not mol.GetAtomWithIdx(x).GetNumRadicalElectrons())
                        and (
                            mol.GetAtomWithIdx(x).GetHybridization()
                            == Chem.HybridizationType.SP3
                        )
                        and (
                            mol.GetAtomWithIdx(x).GetChiralTag()
                            == Chem.ChiralType.CHI_UNSPECIFIED
                        )
                        and (x not in core)
                        for x in (a, b)
                    )
                ):
                    edges.append([a, b])
            sets.append(edges)
            all_edges.extend(edges)
        rec = dict(
            index=g.index,
            members=list(g.atoms),
            heavy_atoms=len(g.heavy_atoms),
            rings=rings,
            eligible_edges=all_edges,
            schemes=[],
            selected_cuts=[],
            status="pending",
        )
        if members & merged:
            rec["status"] = "defer_prior_coarsening"
            records.append(rec)
            continue
        if len(all_edges) > 32:
            rec["status"] = "defer_search_limit"
            records.append(rec)
            continue
        baseline_product, baseline_construction = cap(mol, g.atoms)
        baseline_stereo = stereo_present(baseline_product, baseline_construction)
        accepted = {}
        for edges in sets:
            for pair in combinations(edges, 2):
                ps = split(mol, members, existing, pair)
                sizes = sorted(
                    (
                        sum((mol.GetAtomWithIdx(a).GetAtomicNum() > 1 for a in p))
                        for p in ps
                    )
                )
                reasons = []
                if len(ps) != 2:
                    reasons.append("not_two_components")
                if min(sizes) < 3 or max(sizes) > 20:
                    reasons.append("size_outside3_20")
                if any(
                    (
                        len({i for i, p in enumerate(ps) if set(p) & set(t["atoms"])})
                        > 1
                        for t in motifs
                        if set(t["atoms"]) & members
                    )
                ):
                    reasons.append("core_split")
                allowed = [mol.GetBondBetweenAtoms(*e).GetIdx() for e in pair]
                if not reasons:
                    try:
                        stereo = set()
                        for p in ps:
                            product, construction = cap(
                                mol, p, allowed_ring_bonds=allowed
                            )
                            stereo |= stereo_present(product, construction)
                        if not baseline_stereo <= stereo:
                            reasons.append("new_stereo_loss")
                    except Exception as exc:
                        reasons.append(
                            "cap_failure:" + type(exc).__name__ + ":" + str(exc)
                        )
                scheme = dict(
                    cuts=[list(e) for e in pair],
                    parts=ps,
                    sizes=sizes,
                    rejections=reasons,
                )
                rec["schemes"].append(scheme)
                if not reasons:
                    accepted.setdefault(partition_key(ps), []).append(scheme)
        rec["admissible_partitions"] = len(accepted)
        if len(accepted) == 1:
            options = next(iter(accepted.values()))
            cuts = {tuple(sorted(e)) for s in options for e in s["cuts"]}
            assert len(cuts) == 2
            rec["selected_cuts"] = [list(e) for e in sorted(cuts)]
            rec["status"] = "applied_unique_partition"
            chosen.extend(rec["selected_cuts"])
        else:
            rec["status"] = (
                "defer_no_admissible" if not accepted else "defer_multiple_partitions"
            )
        records.append(rec)
    new = {frozenset(e) for e in chosen}
    finalds = [
        (
            replace(
                d,
                cut=True,
                reason="isolated_macrocycle_unique_pair",
                protections=tuple((p for p in d.protections if p != "ring_bond")),
            )
            if frozenset(d.atoms) in new
            else d
        )
        for d in ds
    ]
    rw = Chem.RWMol(mol)
    for d in finalds:
        if d.cut:
            rw.RemoveBond(*d.atoms)
    gs = [
        ElementaryFragment(
            i,
            tuple(p),
            tuple((a for a in p if mol.GetAtomWithIdx(a).GetAtomicNum() > 1)),
            any((mol.GetAtomWithIdx(a).IsInRing() for a in p)),
        )
        for i, p in enumerate(Chem.GetMolFrags(rw, asMols=False, sanitizeFrags=False))
    ]
    owner = {a: g.index for g in gs for a in g.atoms}
    adj = {g.index: set() for g in gs}
    for d in finalds:
        if d.cut:
            a, b = [owner[i] for i in d.atoms]
            assert a != b
            adj[a].add(b)
            adj[b].add(a)
    return (
        (gs, adj, finalds, motifs, operations),
        dict(blocks=records, new_ring_cuts=chosen),
    )
