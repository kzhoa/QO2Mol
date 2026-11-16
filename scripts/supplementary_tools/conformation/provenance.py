"""Ordered-graph and exact reference-coordinate checks; not canonical identities."""

import hashlib
import json


def _digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def reference_provenance(molecule):
    atoms = [
        dict(
            element=a.elem,
            isotope=a.isotope,
            charge=a.formal_charge,
            bonds=sorted(zip(a.connect, a.bond_type)),
        )
        for a in molecule.Atoms
    ]
    coordinates = [
        [0.0 if float(v) == 0 else float(v) for v in a.coor] for a in molecule.Atoms
    ]
    return dict(
        id=str(molecule.name),
        ordered_graph_sha256=_digest(atoms),
        reference_sha256=_digest(dict(atoms=atoms, coordinates=coordinates)),
    )


def check_scan_reference(molecule, analysis, *, level=None, multiplicity=None):
    if analysis is None:
        return dict(status="unverified", reason="no_scan_analysis")
    name = str(molecule.name)
    report = analysis["report"]
    if name in report.get("failed_molecules", {}):
        raise ValueError(
            f'Scan analysis failed for molecule {name}: {report["failed_molecules"][name]}'
        )
    if name not in analysis["minima"]:
        raise ValueError(f"No successful scan analysis for molecule {name}")
    expected = report.get("provenance", {}).get(name)
    settings = dict(report.get("calculation_settings", {}).get(name, {}))
    if expected:
        settings.update(
            {
                field: expected[field]
                for field in ("level", "multiplicity", "charge")
                if field in expected
            }
        )
    checked = []
    for field, value in (
        ("level", level),
        ("multiplicity", multiplicity),
        ("charge", sum(a.formal_charge for a in molecule.Atoms)),
    ):
        if value is not None and field in settings:
            if settings[field] != value:
                raise ValueError(f"Scan/reference mismatch: {field}")
            checked.append(field)
    if not expected:
        return dict(
            status="unverified",
            reason="historical_scan_without_provenance",
            calculation_settings_checked=checked,
        )
    actual = reference_provenance(molecule)
    for field in ("id", "ordered_graph_sha256", "reference_sha256"):
        if expected.get(field) != actual[field]:
            raise ValueError(f"Scan/reference mismatch: {field}")
    return dict(
        status="verified_reference",
        provenance=expected,
        calculation_settings_checked=checked,
    )
