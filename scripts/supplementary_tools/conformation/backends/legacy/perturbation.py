"""Preserved paired Z-matrix perturbation adapter."""

import copy
import math
import numpy as np
from ...shared.geometry import checked_coordinates, validate_roundtrip
from ...graph import Element
from .source.zmatrix import ZMatrix
from .source.perturb import bond_angle_extend_conformer

from ...shared.constants import ALKANE


def candidate_plan(m, zmat, ignore_alkane):
    """Audit the retained pairing rule; atom indices are zero-based throughout JSON."""

    def same_ring(*ids):
        return any(all(i in r[:-1] for i in ids) for r in m.ring_dict.values())

    def ignored(*ids):
        return ignore_alkane and all(m.Atoms[i].atom_type_name in ALKANE for i in ids)

    bonds = [
        i for i, b in enumerate(zmat.bonds) if b is not None and not ignored(b.a1, b.a2)
    ]
    angles = [
        i
        for i, a in enumerate(zmat.angles)
        if a is not None
        and (a.type21 == "1" or a.type23 == "1")
        and not same_ring(a.a1, a.a3)
        and not ignored(a.a1, a.a2, a.a3)
    ]
    plan = []
    for kind, indices, terms in [("b", bonds, zmat.bonds), ("a", angles, zmat.angles)]:
        half = math.ceil(len(indices) / 2)
        for i, j in zip(indices[:half], (indices + [None])[half : 2 * half]):
            name = f"{kind}-{i+1}" + (f"-{j+1}" if j is not None else "")
            changes = []
            for row, direction in [(i, 1), (j, -1)]:
                if row is None:
                    continue
                term = terms[row]
                if kind == "b":
                    delta = (
                        0.05
                        if same_ring(term.a1, term.a2)
                        else (
                            0.08
                            if term.type == "1"
                            else (0.05 if term.type == "3" else 0.06)
                        )
                    )
                    atoms, reference, unit = [term.a1, term.a2], term.value, "angstrom"
                else:
                    delta, atoms, reference, unit = (
                        10.0,
                        [term.a1, term.a2, term.a3],
                        term.value_a,
                        "degree",
                    )
                changes.append(
                    dict(
                        row=row,
                        atoms=atoms,
                        delta=direction * delta,
                        reference=reference,
                        unit=unit,
                    )
                )
            plan.append(dict(name=name, changes=changes))
    return plan


def generate_q2(m, ignore_alkane):
    m = copy.deepcopy(m)
    for pair in m.Pair14 + m.Pair1n:
        for i in (pair.a1, pair.a2):
            if Element.get(m.Atoms[i].elem).vdw_radius is None:
                raise ValueError(
                    "No original van der Waals radius for " + m.Atoms[i].elem
                )
    reference = np.array([a.coor for a in m.Atoms])
    zmat = ZMatrix()
    zmat.from_molecule(m)
    roundtrip_error = validate_roundtrip(reference, checked_coordinates(zmat))
    plan = candidate_plan(m, zmat, ignore_alkane)
    # Guard the source routine against nonfinite intermediates.
    with np.errstate(all="raise"):
        generated = bond_angle_extend_conformer(m, ignore_alkane=ignore_alkane)
    accepted = {}
    for entry in plan:
        z = generated.get(entry["name"])
        if z is None:
            entry["status"] = "rejected_contact"
            continue
        try:
            accepted[entry["name"]] = checked_coordinates(z)
            entry["status"] = "accepted"
        except (ValueError, FloatingPointError) as exc:
            entry["status"] = "rejected_geometry_adapter"
            entry["error"] = str(exc)
    return accepted, plan, roundtrip_error
