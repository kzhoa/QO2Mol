"""Structured chemical descriptions, distinct from sampling strata."""

from collections import Counter
from rdkit import Chem
from rdkit.Chem import Lipinski
from .rules import recognize


def describe(mol, config):
    heavy = sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
    ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    # Ring-bond connected components include fused, spiro and bridged systems,
    # but exclude non-ring bonds connecting separate systems.
    neighbors = {i: set() for i in ring_atoms}
    for b in mol.GetBonds():
        if b.IsInRing():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            neighbors[i].add(j)
            neighbors[j].add(i)
    remaining = set(ring_atoms)
    systems = []
    while remaining:
        group, frontier = set(), {min(remaining)}
        while frontier:
            group |= frontier
            frontier = {j for i in frontier for j in neighbors[i]} - group
        remaining -= group
        systems.append(len(group))
    elements = Counter(a.GetSymbol() for a in mol.GetAtoms())
    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    motifs = Counter(x["rule"] for x in recognize(mol, config))
    # These deliberately coarse bins do not reproduce legacy_390/legacy_147 classes.
    size_bin = (
        "1-8"
        if heavy <= 8
        else "9-16" if heavy <= 16 else "17-32" if heavy <= 32 else "33+"
    )
    stratum = {
        "size": size_bin,
        "charge": "negative" if charge < 0 else "positive" if charge > 0 else "neutral",
        "rings": "ring" if ring_atoms else "acyclic",
    }
    return {
        "heavy_atoms": heavy,
        "element_counts": dict(sorted(elements.items())),
        "formal_charge": charge,
        "ring_system_sizes": sorted(systems),
        "ring_atoms": len(ring_atoms),
        "aromatic_atoms": sum(a.GetIsAromatic() for a in mol.GetAtoms()),
        "rotatable_bonds": Lipinski.NumRotatableBonds(Chem.RemoveHs(mol)),
        "protected_motif_counts": dict(sorted(motifs.items())),
        "stratum": stratum,
    }
