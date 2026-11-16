"""Graph-local edits, reference-aware clash checking and proper-rotation RMSD."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from ...graph import Element


def positions(molecule):
    xyz = np.asarray([a.coor for a in molecule.Atoms], dtype=float)
    if xyz.shape != (len(molecule.Atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Expected finite N by 3 coordinates")
    return xyz


def unit(vector):
    norm = np.linalg.norm(vector)
    if norm < 1e-10:
        raise ValueError("Coincident atoms or undefined rotation plane")
    return vector / norm


def moving_component(molecule, anchor, endpoint):
    """Atoms on endpoint's side of a bridge, or None for a ring bond."""
    if endpoint not in molecule.Atoms[anchor].connect:
        raise ValueError("The requested edit requires a bonded atom pair")
    seen, pending = {endpoint}, [endpoint]
    while pending:
        i = pending.pop()
        for j in molecule.Atoms[i].connect:
            if {i, j} == {anchor, endpoint}:
                continue
            if j not in seen:
                seen.add(j)
                pending.append(j)
    return None if anchor in seen else sorted(seen)


def angle_degrees(xyz, a, b, c):
    return float(
        np.degrees(
            np.arccos(
                np.clip(np.dot(unit(xyz[a] - xyz[b]), unit(xyz[c] - xyz[b])), -1.0, 1.0)
            )
        )
    )


def rotate(xyz, indices, origin, axis, degrees):
    axis = unit(axis)
    vectors = xyz[indices] - origin
    theta = np.radians(degrees)
    turned = (
        vectors * np.cos(theta)
        + np.cross(axis, vectors) * np.sin(theta)
        + np.outer(vectors @ axis, axis) * (1.0 - np.cos(theta))
    )
    out = xyz.copy()
    out[indices] = turned + origin
    return out


def stretch(xyz, a, b, side, delta):
    if np.linalg.norm(xyz[b] - xyz[a]) + delta <= 1e-8:
        raise ValueError("Non-positive target bond length")
    out = xyz.copy()
    out[side] += delta * unit(xyz[b] - xyz[a])
    return out


def bend(xyz, a, b, c, side, delta):
    target = angle_degrees(xyz, a, b, c) + delta
    if not 1e-5 < target < 180.0 - 1e-5:
        raise ValueError("Target angle must lie strictly between 0 and 180 degrees")
    out = rotate(xyz, side, xyz[b], np.cross(xyz[a] - xyz[b], xyz[c] - xyz[b]), delta)
    if abs(angle_degrees(out, a, b, c) - target) > 1e-6:
        raise ValueError("Angle rotation did not reach its target")
    return out


class ContactCheck:
    """Only reject contacts introduced/worsened relative to the same reference.

    lower_bound = min(scale * (r_i+r_j), reference_distance) - 1e-6 A.
    Pairs separated by at least three graph bonds are examined. Zero disables.
    """

    def __init__(self, molecule, reference, scale):
        pairs = molecule.Pair14 + molecule.Pair1n if scale else []
        self.indices = np.asarray([(p.a1, p.a2) for p in pairs], dtype=int).reshape(
            -1, 2
        )
        radii = []
        for i, j in self.indices:
            ri, rj = (Element.get(molecule.Atoms[k].elem).vdw_radius for k in (i, j))
            if ri is None or rj is None:
                raise ValueError("Missing van der Waals radius for contact check")
            radii.append(ri + rj)
        ref = np.linalg.norm(
            reference[self.indices[:, 0]] - reference[self.indices[:, 1]], axis=1
        )
        self.lower = np.minimum(scale * np.asarray(radii), ref) - 1e-6

    def violation(self, xyz):
        if not np.isfinite(xyz).all():
            raise ValueError("Non-finite candidate geometry")
        distances = np.linalg.norm(
            xyz[self.indices[:, 0]] - xyz[self.indices[:, 1]], axis=1
        )
        bad = np.flatnonzero(distances < self.lower)
        if not len(bad):
            return None
        i = int(bad[0])
        return dict(
            atoms=self.indices[i].tolist(),
            distance=float(distances[i]),
            lower_bound=float(self.lower[i]),
        )


def aligned_rmsd(first, second, indices):
    """Fixed atom correspondence, translation/proper rotation only; no reflection."""
    p, q = np.asarray(first)[indices], np.asarray(second)[indices]
    if len(p) == 0:
        raise ValueError("No atoms selected for RMSD")
    p, q = p - p.mean(axis=0), q - q.mean(axis=0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.diag([1.0, 1.0, 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0])
    return float(np.sqrt(np.mean(np.sum((p @ (u @ correction @ vt) - q) ** 2, axis=1))))


def rdkit_reference(molecule, xyz):
    rd = molecule.to_rdkit()
    rd.RemoveAllConformers()
    conformer = Chem.Conformer(len(molecule.Atoms))
    conformer.Set3D(True)
    for i, point in enumerate(xyz):
        conformer.SetAtomPosition(i, point.tolist())
    rd.AddConformer(conformer)
    return rd


def set_torsions(reference_rd, terms, angles):
    """Set each acyclic torsion on a fresh conformer and validate final targets."""
    rd = Chem.Mol(reference_rd)
    conformer = rd.GetConformer()
    for atoms, angle in zip(terms, angles):
        rdMolTransforms.SetDihedralDeg(conformer, *atoms, float(angle))
    xyz = np.asarray(conformer.GetPositions())
    if not np.isfinite(xyz).all():
        raise ValueError("Non-finite dihedral rotation")
    for atoms, angle in zip(terms, angles):
        achieved = rdMolTransforms.GetDihedralDeg(conformer, *atoms)
        if (
            not np.isfinite(achieved)
            or abs((achieved - angle + 180.0) % 360.0 - 180.0) > 1e-5
        ):
            raise ValueError("Combination did not preserve all requested torsions")
    return xyz
