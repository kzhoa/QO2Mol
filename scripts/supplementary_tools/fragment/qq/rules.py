"""Overlapping structural recognition, separate from explicit bond decisions."""

from functools import lru_cache
from rdkit import Chem
from .models import BondDecision

# QQ-specific design rules, not a reconstruction of the legacy_147/legacy_390 tables.
PROTECTED = (
    ("acyl_heteroatom", "[C,S,P](=[O,S,N])-[N,O,S]"),
    ("amidinium", "[C](=[N])-[N]"),
    ("nitro", "[N+](=[O])-[O-]"),
)


@lru_cache(maxsize=128)
def query(smarts):
    result = Chem.MolFromSmarts(smarts)
    if result is None:
        raise ValueError("Invalid protection SMARTS: " + smarts)
    return result


def recognize(mol, config):
    matches = []
    rules = PROTECTED + tuple(
        (f"custom_{i}", s) for i, s in enumerate(config.protected_smarts)
    )
    for name, smarts in rules:
        pattern = query(smarts)
        # Fail explicitly rather than silently truncate motif matches at RDKit's
        # usual default limit; the bound applies independently to each pattern.
        found = mol.GetSubstructMatches(pattern, uniquify=True, maxMatches=100001)
        if len(found) > 100000:
            raise ValueError("Protection match limit exceeded: " + name)
        for atoms in found:
            bonds = tuple(
                mol.GetBondBetweenAtoms(
                    atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()]
                ).GetIdx()
                for b in pattern.GetBonds()
            )
            matches.append({"rule": name, "atoms": list(atoms), "bonds": list(bonds)})
    return matches


def decide_bonds(mol, config):
    motifs = recognize(mol, config)
    protected = {}
    for motif in motifs:
        for b in motif["bonds"]:
            protected.setdefault(b, set()).add(motif["rule"])
    decisions = []
    degree = {
        a.GetIdx(): sum(n.GetAtomicNum() > 1 for n in a.GetNeighbors())
        for a in mol.GetAtoms()
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        ring_link = a.IsInRing() and b.IsInRing() and not bond.IsInRing()
        protections = set(protected.get(bond.GetIdx(), ()))
        if min(a.GetAtomicNum(), b.GetAtomicNum()) == 1:
            protections.add("hydrogen_attachment")
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
            protections.add("not_nonaromatic_single")
        if bond.IsInRing():
            protections.add("ring_bond")
        if config.protect_charged and (a.GetFormalCharge() or b.GetFormalCharge()):
            protections.add("charged_endpoint")
        if config.protect_terminal and min(degree[a.GetIdx()], degree[b.GetIdx()]) <= 1:
            protections.add("terminal_group")
        if config.protect_conjugated and bond.GetIsConjugated() and not ring_link:
            protections.add("conjugated_bond")
        candidates = []
        if ring_link:
            candidates.append("inter_ring_link")
        if a.IsInRing() != b.IsInRing():
            candidates.append("ring_sidechain_link")
        if a.GetAtomicNum() not in (1, 6) or b.GetAtomicNum() not in (1, 6):
            candidates.append("heteroatom_link")
        if max(degree[a.GetIdx()], degree[b.GetIdx()]) >= 3:
            candidates.append("branch_link")
        if config.cut_acyclic_cc and a.GetAtomicNum() == b.GetAtomicNum() == 6:
            candidates.append("acyclic_cc")
        cut = bool(candidates) and not protections
        reason = (
            "protected" if protections else ("eligible_link" if cut else "no_cut_rule")
        )
        decisions.append(
            BondDecision(
                bond.GetIdx(),
                (a.GetIdx(), b.GetIdx()),
                cut,
                reason,
                tuple(sorted(protections)),
                tuple(candidates),
            )
        )
    return decisions, motifs
