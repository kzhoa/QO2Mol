"""Legacy conformational sampling backend."""

from .torsion import generate_q1
from .q4 import generate_q4, minima_from_scan_records
from .perturbation import generate_q2 as _generate_q2


def generate_q2(molecule, ignore_alkane):
    coordinates, plan, error = _generate_q2(molecule, ignore_alkane)
    return coordinates, dict(candidate_plan=plan, roundtrip_max_error_angstrom=error)
