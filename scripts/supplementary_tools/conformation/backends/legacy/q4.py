"""Optional Q4 combinations, explicit optimization hook and legacy deduplication."""

import copy
import itertools
import math
import types
import numpy as np

from .q4_geometry import rotating_side, MP
from .source.combination import (
    create_lm_by_combine_scan_rlm,
    remove_similar_conformer,
    local_minimum_pes,
)


def _update_topol_value(self):
    # The copied generator calls this before rotation, but its returned topology
    # values are never consumed; the rotation adapter reads the live coordinates.
    pass


def minima_from_scan_records(records):
    """Use saved conform_type labels, sort minima by energy; never invent labels.

    Input: id, torsion (four zero-based indices), angle (degrees), energy,
    conform_type. Energies must have consistent units and calculation settings.
    Returns id -> 'i-j-k-l' -> angles. No coordinate/reference selection here.
    """
    curves = {}
    for r in records:
        if not math.isfinite(float(r["energy"])) or not math.isfinite(
            float(r["angle"])
        ):
            raise ValueError("Finite scan energies/angles required")
        if (
            len(r["torsion"]) != 4
            or len(set(r["torsion"])) != 4
            or any(type(i) is not int or i < 0 for i in r["torsion"])
        ):
            raise ValueError("Four zero-based torsion indices required")
        term = "-".join(map(str, r["torsion"]))
        curves.setdefault(str(r["id"]), {}).setdefault(term, []).append(
            dict(
                energy=float(r["energy"]),
                conform_type=r["conform_type"],
                constrain=[[r["torsion"], float(r["angle"])]],
            )
        )
    return local_minimum_pes(curves) if curves else {}


def real_values(values):
    values = np.asarray(values)
    if not np.isfinite(values).all() or np.max(np.abs(values.imag)) > 1e-8:
        raise ValueError("Non-real/nonfinite inertia eigenvalues")
    return values.real.tolist()


def generate_q4(molecule, minima, *, limit=64, optimizer=None):
    """Return retained coordinates and full audit; optimizer(m) may mutate m.

    The CLI chooses optimizer=None explicitly, matching optimize_structure=False.
    A supplied callback must preserve atom order, graph, charge and atom identity.
    """
    if type(limit) is not int or limit < 1:
        raise ValueError("Q4 limit must be a positive integer")
    if not isinstance(minima, dict) or not minima:
        raise ValueError("Q4 requires an ordered nonempty torsion-to-minima mapping")
    reference = copy.deepcopy(molecule)
    for term, values in minima.items():
        torsion = [int(s) for s in term.split("-")]
        rotating_side(reference, torsion)
        if not isinstance(values, list) or any(
            not isinstance(v, (int, float)) or not math.isfinite(v) for v in values
        ):
            raise ValueError("Q4 minima must be finite angle lists in degrees")
    reference.update_topol_value = types.MethodType(_update_topol_value, reference)
    with np.errstate(all="raise"):
        candidates = create_lm_by_combine_scan_rlm(
            reference, minima, n=limit, create_constrain=False
        )
    changes = list(itertools.islice(itertools.product(*minima.values()), limit))
    for m in candidates:
        if optimizer is not None:
            original_graph = [
                (a.elem, a.formal_charge, list(a.connect), list(a.bond_type))
                for a in m.Atoms
            ]
            optimizer(m)
            if [
                (a.elem, a.formal_charge, list(a.connect), list(a.bond_type))
                for a in m.Atoms
            ] != original_graph:
                raise ValueError("Q4 optimizer changed atom order or molecular graph")
        if not np.isfinite(np.asarray([a.coor for a in m.Atoms])).all():
            raise ValueError("Q4 optimizer/rotation produced nonfinite coordinates")
    retained = remove_similar_conformer(candidates, target_molecule=reference)
    keep_ids = {id(m) for m in retained}
    reference_inertia = MP._inertia_calculate(reference)[0]
    plan, coordinates = [], {}
    for i, (m, angles) in enumerate(zip(candidates, changes)):
        name = f"combination-{i+1}"
        inertia = MP._inertia_calculate(m)[0]
        kept = id(m) in keep_ids
        plan.append(
            dict(
                name=name,
                angles=dict(zip(minima, angles)),
                status="accepted" if kept else "rejected_legacy_inertia",
                inertia=real_values(inertia),
            )
        )
        if kept:
            coordinates[name] = np.asarray([a.coor for a in m.Atoms])
    return coordinates, dict(
        candidate_plan=plan,
        total_combinations=math.prod(len(v) for v in minima.values()),
        generated_count=len(candidates),
        accepted_count=len(retained),
        limit=limit,
        reference_inertia=real_values(reference_inertia),
        optimization="none" if optimizer is None else "external_callback",
        deduplication="legacy wrapper: includes H; uncentered, unsorted inertia eigenvalues; atol=5, rtol=0",
    )
