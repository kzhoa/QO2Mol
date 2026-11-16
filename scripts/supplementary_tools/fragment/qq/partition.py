"""Connected components after the declared QQ cuts, with explicit adjacency."""

from rdkit import Chem
from .models import ElementaryFragment
from .rules import decide_bonds


def _baseline_partition(mol, config):
    decisions, motifs = decide_bonds(mol, config)
    graph = Chem.RWMol(mol)
    for decision in decisions:
        if decision.cut:
            graph.RemoveBond(*decision.atoms)
    groups = Chem.GetMolFrags(graph, asMols=False, sanitizeFrags=False)
    fragments = [
        ElementaryFragment(
            i,
            tuple(g),
            tuple(a for a in g if mol.GetAtomWithIdx(a).GetAtomicNum() > 1),
            any(mol.GetAtomWithIdx(a).IsInRing() for a in g),
        )
        for i, g in enumerate(groups)
    ]
    owners = {a: f.index for f in fragments for a in f.atoms}
    adjacency = {f.index: set() for f in fragments}
    for decision in decisions:
        if decision.cut:
            i, j = map(owners.get, decision.atoms)
            if i == j:
                raise ValueError("Declared cut did not separate components")
            adjacency[i].add(j)
            adjacency[j].add(i)
    return fragments, adjacency, decisions, motifs


def partition(mol, config):
    """Public partition is the same basic-units algorithm as process()."""
    from .basic_units import partition_basic_units
    return partition_basic_units(mol, config)[:4]
