"""Stable document-based interface with explicit legacy/qq backend selection."""

from importlib import import_module
from .inputs import from_document
from .provenance import check_scan_reference
from .preparation import prepare_geometry
from .gaussian import gaussian_input, molecular_charge, validate_level
from .options import Q2Options, Q4Options, ScanOptions

BACKENDS = ("legacy", "qq")


def get_backend(backend):
    if backend not in BACKENDS:
        raise ValueError("backend must be legacy or qq")
    return import_module(f".backends.{backend}", __package__)


def resolve_options(backend, stage, options):
    get_backend(backend)
    if backend == "legacy":
        if options is not None:
            raise ValueError("qq options cannot be applied to the legacy backend")
        return {}
    cls = {"q2": Q2Options, "q4": Q4Options}[stage]
    if options is not None and not isinstance(options, cls):
        raise TypeError(f"{stage} requires {cls.__name__}")
    return {"options": options or cls()}


def scan_torsions(
    document, *, alkane_policy="supplied-types", backend="qq", return_report=False
):
    engine = get_backend(backend)
    m = prepare_geometry(
        from_document(document),
        document.get("atom_types"),
        alkane_policy,
        backend=backend,
        stage="q1",
    )
    if return_report:
        if backend != "qq":
            raise ValueError("Q1 decision reports require qq")
        return engine.generate_q1(
            m, alkane_policy == "supplied-types", return_report=True
        )
    return engine.generate_q1(m, alkane_policy == "supplied-types")


def perturb_geometry(
    document, *, alkane_policy="supplied-types", backend="qq", options=None
):
    engine = get_backend(backend)
    arguments = resolve_options(backend, "q2", options)
    m = prepare_geometry(
        from_document(document),
        document.get("atom_types"),
        alkane_policy,
        backend=backend,
        stage="q2",
    )
    coordinates, report = engine.generate_q2(
        m, alkane_policy == "supplied-types", **arguments
    )
    return dict(
        coordinates={k: v.tolist() for k, v in coordinates.items()},
        backend=backend,
        **report,
    )


def combine_minima(
    document,
    minima,
    *,
    limit=64,
    optimizer=None,
    backend="qq",
    options=None,
    scan_analysis=None,
):
    engine = get_backend(backend)
    arguments = resolve_options(backend, "q4", options)
    m = prepare_geometry(
        from_document(document), alkane_policy="keep", backend=backend, stage="q4"
    )
    if scan_analysis is not None and backend != "qq":
        raise ValueError("Scan provenance checking requires qq")
    verification = check_scan_reference(m, scan_analysis) if backend == "qq" else None
    if scan_analysis is not None and minima != scan_analysis["minima"][str(m.name)]:
        raise ValueError("Minima differ from supplied scan analysis")
    coordinates, report = engine.generate_q4(
        m, minima, limit=limit, optimizer=optimizer, **arguments
    )
    if verification is not None:
        report["scan_reference_check"] = verification
    return dict(
        coordinates={k: v.tolist() for k, v in coordinates.items()},
        backend=backend,
        **report,
    )


def optimization_input(document, *, level, multiplicity, backend="qq"):
    get_backend(backend)
    validate_level(level)
    m = prepare_geometry(
        from_document(document), alkane_policy="keep", backend=backend, stage="q0"
    )
    charge = molecular_charge(m, multiplicity)
    return gaussian_input(
        m,
        [a.coor for a in m.Atoms],
        charge,
        multiplicity,
        level,
        "Q0",
        "reference",
        backend=backend,
    )


def minima_from_scan_records(
    records, *, backend="qq", options=None, return_report=False
):
    engine = get_backend(backend)
    if backend == "legacy":
        if options is not None or return_report:
            raise ValueError("Scan options/reports require qq")
        return engine.minima_from_scan_records(records)
    if options is not None and not isinstance(options, ScanOptions):
        raise TypeError("ScanOptions required")
    return engine.minima_from_scan_records(
        records, options=options, return_report=return_report
    )
