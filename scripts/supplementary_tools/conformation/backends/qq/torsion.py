"""Historical flexible-bond eligibility with geometry-aware scan definitions."""

import copy
import numpy as np
from rdkit import Chem
from ...shared.torsion import flexible_torsion
from ...shared.constants import ALKANE
from ...provenance import reference_provenance
from .geometry import moving_component, positions


def generate_q1(molecule, ignore_alkane, *, return_report=False):
    m = copy.deepcopy(molecule)
    for bond in m.Bonds:
        bond.flexible = "no"
    flexible_torsion(m)
    xyz = positions(m)
    ranks = list(Chem.CanonicalRankAtoms(m.to_rdkit(), breakTies=False))
    torsions, decisions = [], []
    for bond in sorted(m.Bonds, key=lambda b: (b.a1, b.a2)):
        b, c = bond.a1, bond.a2
        entry = dict(bond=[b, c])
        decisions.append(entry)
        if getattr(bond, "flexible", "no") != "yes":
            entry["reason"] = "historical_ineligible"
            continue
        if moving_component(m, b, c) is None:
            entry["reason"] = "ring_bond"
            continue
        if any("3" in m.Atoms[i].bond_type for i in (b, c)):
            entry["reason"] = "adjacent_triple_bond"
            continue
        axis = xyz[c] - xyz[b]
        norm = np.linalg.norm(axis)
        if norm < 1e-8:
            entry["reason"] = "undefined_axis"
            continue
        axis /= norm
        chosen = []
        for center, other in ((b, c), (c, b)):
            candidates = []
            for i in m.Atoms[center].connect:
                if i == other or m.Atoms[i].s == "dummy":
                    continue
                vector = xyz[i] - xyz[center]
                length = np.linalg.norm(vector)
                quality = (
                    np.linalg.norm(np.cross(vector / length, axis))
                    if length > 1e-8
                    else 0.0
                )
                if quality > 1e-4:
                    # Prefer heavy references, then well-conditioned planes;
                    # topology rank and index break exact geometric ties.
                    candidates.append(
                        (
                            m.Atoms[i].elem != "H",
                            round(float(quality), 10),
                            ranks[i],
                            -i,
                            i,
                        )
                    )
            if candidates:
                chosen.append(max(candidates)[-1])
        if len(chosen) != 2:
            entry["reason"] = "undefined_dihedral_plane"
            continue
        term = [chosen[0], b, c, chosen[1]]
        entry["torsion"] = term
        if ignore_alkane and all(m.Atoms[i].atom_type_name in ALKANE for i in term):
            entry["reason"] = "alkane_filter"
            continue
        entry["reason"] = "accepted"
        torsions.append(term)
    return (
        dict(
            torsions=torsions,
            bond_decisions=decisions,
            provenance=dict(reference_provenance(m), expected_torsions=torsions),
        )
        if return_report
        else torsions
    )
