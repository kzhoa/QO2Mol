"""Shared input validation and topology adaptation for both backends."""

import copy
from types import SimpleNamespace
import numpy as np
from rdkit import Chem
from .graph import prepare
from .shared.geometry import Bond
from .graph import Element


def prepare_geometry(
    molecule,
    atom_types=None,
    alkane_policy="supplied-types",
    *,
    backend="legacy",
    stage=None
):
    if backend not in {"legacy", "qq"} or stage not in {None, "q0", "q1", "q2", "q4"}:
        raise ValueError("Unknown preparation backend/stage")
    m = copy.deepcopy(molecule)
    coords = np.asarray([a.coor for a in m.Atoms], dtype=float)
    if (
        len(m.Atoms) < 2
        or coords.shape != (len(m.Atoms), 3)
        or not np.isfinite(coords).all()
    ):
        raise ValueError(
            "At least two atoms with finite explicit coordinates are required"
        )
    rd = m.to_rdkit()
    if len(Chem.GetMolFrags(rd)) != 1:
        raise ValueError("One connected molecule per record is required")
    if any(a.GetNumImplicitHs() for a in rd.GetAtoms()):
        raise ValueError("Input must contain explicit hydrogen atoms")
    # to_rdkit suppresses implicit H: also check valence using a normal sanitized copy.
    check = Chem.Mol(rd)
    for a in check.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(check)
    if any(a.GetNumImplicitHs() for a in check.GetAtoms()):
        raise ValueError("Input must contain explicit hydrogen atoms")
    if alkane_policy not in {"keep", "supplied-types"}:
        raise ValueError("Unknown alkane policy")
    if atom_types is not None and (
        len(atom_types) != len(m.Atoms)
        or not all(isinstance(t, str) and t for t in atom_types)
    ):
        raise ValueError(
            "atom_types must contain one nonempty historical type per atom"
        )
    if alkane_policy == "supplied-types" and atom_types is None:
        raise ValueError(
            "Historical atom_types required; explicitly use --alkane-policy keep if unavailable"
        )
    if backend == "legacy" or stage in {None, "q1"}:
        prepare(m)
    else:
        m.create_topols()
    for i, atom in enumerate(m.Atoms):
        atom.s = "atom"
        atom.mass = Element.get(atom.elem).mass
        atom.atom_type_name = atom_types[i] if atom_types is not None else ""
    m.Bonds = [Bond(b.a1, b.a2) for b in m.Bonds]
    for b in m.Bonds:
        b.calc_value([m.Atoms[b.a1].coor, m.Atoms[b.a2].coor])
    if backend == "qq" and stage in {"q0", "q1"}:
        return m
    distances = Chem.GetDistanceMatrix(rd)
    m.Pair14, m.Pair1n = [], []
    for i in range(len(m.Atoms)):
        for j in range(i + 1, len(m.Atoms)):
            if distances[i, j] >= 3:
                group = m.Pair14 if distances[i, j] == 3 else m.Pair1n
                group.append(SimpleNamespace(a1=i, a2=j))
    return m
