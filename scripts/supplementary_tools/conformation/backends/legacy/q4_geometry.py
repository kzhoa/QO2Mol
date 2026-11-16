"""Topology adapters for retained Q4 combination and legacy inertia behavior."""

import numpy as np
from .source.rotation import calculate_dihedral, rotation_dihedral
from .source.inertia import calc_inertia


def rotating_side(m, torsion):
    if (
        len(torsion) != 4
        or len(set(torsion)) != 4
        or any(i < 0 or i >= len(m.Atoms) for i in torsion)
    ):
        raise ValueError("Q4 torsion must contain four distinct valid atom indices")
    if any(b not in m.Atoms[a].connect for a, b in zip(torsion, torsion[1:])):
        raise ValueError("Q4 torsion atoms must form a bonded path")
    _, b, c, _ = torsion
    seen, pending = {c}, [c]
    while pending:
        i = pending.pop()
        for j in m.Atoms[i].connect:
            if {i, j} == {b, c}:
                continue
            if j not in seen:
                seen.add(j)
                pending.append(j)
    if b in seen:
        raise ValueError(
            "Q4 central bond is in a ring; rigid side rotation is unsupported"
        )
    return sorted(seen - {c})


class ME:
    @staticmethod
    def _structure_change(molecule, atoms, value):
        side = rotating_side(molecule, atoms)
        coords = [molecule.Atoms[i].coor for i in atoms]
        with np.errstate(all="raise"):
            delta = value - calculate_dihedral(*coords)
            for i in side:
                molecule.Atoms[i].coor = rotation_dihedral(
                    coords[1], coords[2], molecule.Atoms[i].coor, delta
                )
            result = calculate_dihedral(*[molecule.Atoms[i].coor for i in atoms])
        if not np.isfinite(result) or abs((result - value + 180) % 360 - 180) > 1e-5:
            raise ValueError("Q4 rotation failed to reach requested dihedral")


class MP:
    @staticmethod
    def _inertia_calculate(molecule, ignore_hydrogen=False):
        # Retain actual source wrapper behavior, including its ignored argument.
        return calc_inertia(molecule, ignore_hydrogen=False)
