"""New compatibility adapter: ZMatrix angles in radians, value_a in degrees."""

import numpy as np


def unit(v):
    length = np.linalg.norm(v)
    if length < 1e-10:
        raise ValueError("Coincident or collinear Z-matrix references")
    return v / length


class Bond:
    def __init__(self, a1, a2):
        self.a1, self.a2 = a1, a2

    def calc_value(self, coords):
        self.value = float(np.linalg.norm(np.asarray(coords[0]) - coords[1]))
        if self.value < 1e-10:
            raise ValueError("Zero bond length")


class Angle:
    def __init__(self, style, a1, a2, a3):
        self.a1, self.a2, self.a3 = a1, a2, a3

    def calc_value(self, coords):
        a, b, c = np.asarray(coords, dtype=float)
        self.value = float(np.arccos(np.clip(np.dot(unit(a - b), unit(c - b)), -1, 1)))
        self.value_a = float(np.degrees(self.value))


class Dihedral:
    def __init__(self, style, a1, a2, a3, a4):
        self.a1, self.a2, self.a3, self.a4 = a1, a2, a3, a4

    def calc_value(self, coords):
        # Inverse of the retained ZMatrix.to_coordinates() local basis.
        d, c, b, a = np.asarray(coords, dtype=float)
        bc = c - b
        n = unit(np.cross(d - c, bc))
        m = unit(np.cross(bc, n))
        self.value = float(np.arctan2(np.dot(a - b, n), np.dot(a - b, m)))
        self.value_a = float(np.degrees(self.value))


def checked_coordinates(zmat):
    # Earlier references and non-linear angles are required by the legacy basis.
    for i, angle in enumerate(zmat.angles[: len(zmat.atoms)]):
        if angle is not None and not 1e-7 < angle.value < np.pi - 1e-7:
            raise ValueError("Linear or out-of-range Z-matrix angle")
        if angle is not None and not (angle.a1 < i and angle.a2 < i):
            raise ValueError("Forward Z-matrix angle reference")
    for i, dih in enumerate(zmat.dihedrals[: len(zmat.atoms)]):
        if dih is not None and not all(j < i for j in (dih.a1, dih.a2, dih.a3)):
            raise ValueError("Forward Z-matrix dihedral reference")
    with np.errstate(all="raise"):
        coords = zmat.to_coordinates()
    if not np.isfinite(coords).all():
        raise ValueError("Non-finite reconstructed coordinates")
    return coords


def validate_roundtrip(original, reconstructed):
    original, reconstructed = np.asarray(original), np.asarray(reconstructed)
    # Proper-rotation Kabsch alignment also rejects mirror inversion.
    p = original - original.mean(axis=0)
    q = reconstructed - reconstructed.mean(axis=0)
    u, _, vt = np.linalg.svd(p.T @ q)
    sign = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1.0, 1.0, sign]) @ vt
    error = float(np.max(np.linalg.norm(p @ rotation - q, axis=1)))
    if error > 1e-6:
        raise ValueError(
            f"Z-matrix roundtrip changed geometry/chirality: {error:.6g} A"
        )
    return error
