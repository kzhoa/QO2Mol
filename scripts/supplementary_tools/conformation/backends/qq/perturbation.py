"""Independent symmetric bond/angle edits, without an atom-order Z-matrix."""

from dataclasses import asdict
from itertools import combinations
import numpy as np
from ...options import Q2Options
from ...shared.constants import ALKANE
from .geometry import (
    positions,
    moving_component,
    stretch,
    bend,
    angle_degrees,
    ContactCheck,
)


def generate_q2(molecule, ignore_alkane, *, options=None):
    settings = options or Q2Options()
    reference = positions(molecule)
    contact = ContactCheck(molecule, reference, settings.clash_scale)
    candidates, plan, skipped = {}, [], []

    def ignored(atoms):
        return ignore_alkane and all(
            molecule.Atoms[i].atom_type_name in ALKANE for i in atoms
        )

    def emit(kind, atoms, side, value, step, operation):
        for sign, suffix in [(-1, "minus"), (1, "plus")]:
            name = kind + "-" + "-".join(map(str, atoms)) + "-" + suffix
            entry = dict(
                name=name,
                atoms=list(atoms),
                moving_atoms=side,
                reference=value,
                delta=sign * step,
                target=value + sign * step,
                unit="angstrom" if kind == "bond" else "degree",
            )
            plan.append(entry)
            try:
                xyz = operation(sign * step)
                violation = contact.violation(xyz)
                if violation is not None:
                    entry.update(status="rejected_contact", contact=violation)
                else:
                    candidates[name] = xyz
                    entry["status"] = "accepted"
            except (ValueError, FloatingPointError) as exc:
                entry.update(status="rejected_geometry", error=str(exc))

    for bond in sorted(molecule.Bonds, key=lambda b: (b.a1, b.a2)):
        a, b = bond.a1, bond.a2
        if ignored([a, b]):
            skipped.append(dict(kind="bond", atoms=[a, b], reason="alkane_filter"))
            continue
        side = moving_component(molecule, a, b)
        if side is None:
            skipped.append(dict(kind="bond", atoms=[a, b], reason="ring_closure"))
            continue
        value = float(np.linalg.norm(reference[b] - reference[a]))
        emit(
            "bond",
            [a, b],
            side,
            value,
            settings.bond_step,
            lambda delta: stretch(reference, a, b, side, delta),
        )

    for b, atom in enumerate(molecule.Atoms):
        for a, c in combinations(sorted(atom.connect), 2):
            atoms = [a, b, c]
            if ignored(atoms):
                skipped.append(dict(kind="angle", atoms=atoms, reason="alkane_filter"))
                continue
            if not any(atom.bond_type[atom.connect.index(i)] == "1" for i in (a, c)):
                skipped.append(
                    dict(kind="angle", atoms=atoms, reason="no_adjacent_single_bond")
                )
                continue
            side = moving_component(molecule, b, c)
            if side is None:
                side = moving_component(molecule, b, a)
                if side is None:
                    skipped.append(
                        dict(kind="angle", atoms=atoms, reason="ring_closure")
                    )
                    continue
                # Swap endpoints so bend() always rotates its third-atom side.
                i, j = c, a
            else:
                i, j = a, c
            try:
                value = angle_degrees(reference, i, b, j)
            except ValueError as exc:
                skipped.append(
                    dict(
                        kind="angle",
                        atoms=atoms,
                        reason="undefined_geometry",
                        error=str(exc),
                    )
                )
                continue
            emit(
                "angle",
                [i, b, j],
                side,
                value,
                settings.angle_step,
                lambda delta: bend(reference, i, b, j, side, delta),
            )
    return candidates, dict(
        candidate_plan=plan,
        skipped_coordinates=skipped,
        proposed_count=len(plan),
        accepted_count=len(candidates),
        settings=asdict(settings),
        scheme="independent minus/plus edits from a common reference; bridge-side motion",
    )
