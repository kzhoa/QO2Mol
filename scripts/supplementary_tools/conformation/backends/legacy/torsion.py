"""Preserved flexible torsion selection."""

import copy
from ...shared.torsion import flexible_torsion, scan_torsion
from .source.perturb import ignore_alkane_torsion


def generate_q1(m, ignore_alkane):
    m = copy.deepcopy(m)
    flexible_torsion(m)
    scan_torsion(m)
    if ignore_alkane:
        ignore_alkane_torsion(m)
    return m.torsions
