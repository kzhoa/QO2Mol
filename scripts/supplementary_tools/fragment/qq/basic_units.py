"""Validated one-pass bounded coarsening stage of qq_basic_units_v1.

The initial bond partition is an internal stage, not a selectable backend.
"""
from dataclasses import replace
from rdkit import Chem
from .models import ElementaryFragment
from .partition import _baseline_partition as partition


def partition_basic_units(mol, config):
    """Return fragments, adjacency, decisions, motifs, and auditable proposals.

Callers must use read_parent first, as for partition(). The basic-units contract is
limited to the supported elementary configuration; incompatible partition options
are rejected rather than silently ignored. All proposals use original blocks.
    """
    if (config.ether_linker_max_heavy or config.cut_acyclic_cc
            or not config.protect_conjugated or not config.protect_charged
            or not config.protect_terminal):
        raise ValueError('Basic-unit merging requires the supported boundary and protection options')
    groups, _, decisions, motifs = partition(mol, config)
    owner = {a: g.index for g in groups for a in g.atoms}
    core = {a for motif in motifs for a in motif.get('atoms', ())}
    proposals = []
    for atom in mol.GetAtoms():
        center = atom.GetIdx()
        kind = {8: 'ether', 7: 'amine', 6: 'branch'}.get(atom.GetAtomicNum())
        neighbors = [a for a in atom.GetNeighbors() if a.GetAtomicNum() > 1]
        degree = len(neighbors)
        if (not kind or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsAromatic()
                or center in core or atom.GetHybridization() != Chem.HybridizationType.SP3
                or not ((kind == 'ether' and degree == 2)
                        or (kind == 'amine' and degree in (2, 3))
                        or (kind == 'branch' and degree in (3, 4)))):
            continue
        if any(a.GetAtomicNum() != 6 or a.IsInRing() or a.GetFormalCharge()
               or a.GetHybridization() != Chem.HybridizationType.SP3 or a.GetIdx() in core
               for a in neighbors):
            continue
        shell = {center, *(a.GetIdx() for a in neighbors)}
        support = {owner[a] for a in shell}
        if len(support) == 1:
            continue
        members = {a for g in groups if g.index in support for a in g.atoms}
        size = sum(mol.GetAtomWithIdx(a).GetAtomicNum() > 1 for a in members)
        proposals.append(dict(kind=kind, center=center, shell=sorted(shell),
                              baseline_blocks=sorted(support), members=sorted(members),
                              heavy_atoms=size))
    for i, p in enumerate(proposals):
        p['conflicts'] = [j for j, q in enumerate(proposals) if i != j
                          and set(p['baseline_blocks']) & set(q['baseline_blocks'])]
        p['status'] = ('defer_shared_baseline_block' if p['conflicts'] else
                       'defer_large_union' if p['heavy_atoms'] > 20 else 'applied')
        p['cancelled_bonds'] = []
    merged_owner = dict(owner)
    selected = {}
    for i, p in enumerate(proposals):
        if p['status'] == 'applied':
            for a in p['members']:
                merged_owner[a] = min(p['baseline_blocks'])
                selected[a] = i
    final_decisions = []
    for d in decisions:
        a, b = d.atoms
        if d.cut and merged_owner[a] == merged_owner[b]:
            p = proposals[selected[a]]
            p['cancelled_bonds'].append(list(d.atoms))
            reason = 'basic_unit_' + p['kind']
            d = replace(d, cut=False, reason=reason, protections=(*d.protections, reason))
        final_decisions.append(d)
    # Build final connected components from the revised cut list, independently
    # of the proposed ownership map; this also verifies complete block unions.
    graph = Chem.RWMol(mol)
    for d in final_decisions:
        if d.cut:
            graph.RemoveBond(*d.atoms)
    components = Chem.GetMolFrags(graph, asMols=False, sanitizeFrags=False)
    expected = {}
    for a, label in merged_owner.items():
        expected.setdefault(label, set()).add(a)
    if {frozenset(g) for g in components} != {frozenset(g) for g in expected.values()}:
        raise ValueError('Cancelled cuts do not reproduce complete proposed block unions')
    fragments = [ElementaryFragment(i, tuple(g),
                 tuple(a for a in g if mol.GetAtomWithIdx(a).GetAtomicNum() > 1),
                 any(mol.GetAtomWithIdx(a).IsInRing() for a in g))
                 for i, g in enumerate(components)]
    final_owner = {a: g.index for g in fragments for a in g.atoms}
    adjacency = {g.index: set() for g in fragments}
    for d in final_decisions:
        if d.cut:
            a, b = (final_owner[x] for x in d.atoms)
            if a == b:
                raise ValueError('Declared cut did not separate components')
            adjacency[a].add(b)
            adjacency[b].add(a)
    return fragments, adjacency, final_decisions, motifs, proposals
