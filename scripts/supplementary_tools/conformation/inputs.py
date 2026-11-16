"""Validate conformation documents before the historical graph adapter coerces values."""

from numbers import Integral, Real
import numpy as np
from .graph import from_document as graph_from_document


def from_document(document):
    n = len(document["elements"])
    for field in ("connectivity", "bond_type", "formal_charge", "coordinates"):
        if field not in document or len(document[field]) != n:
            raise ValueError(f"{field} must contain exactly one entry per atom")
    for value in document["formal_charge"]:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise ValueError("formal_charge values must be integers, without coercion")
    for row in document["connectivity"]:
        if any(isinstance(i, bool) or not isinstance(i, Integral) for i in row):
            raise ValueError("Connectivity indices must be integers")
    for row in document["coordinates"]:
        if len(row) != 3 or any(
            isinstance(v, bool) or not isinstance(v, Real) for v in row
        ):
            raise ValueError("Coordinates must be numeric N x 3 values")
    if not np.isfinite(np.asarray(document["coordinates"], dtype=float)).all():
        raise ValueError("Coordinates must be finite")
    return graph_from_document(document)
