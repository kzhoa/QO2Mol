"""Prepare Q0/Q1/Q2/Q4 tasks using the selected conformation backend."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
from dataclasses import asdict
from pathlib import Path
import numpy as np
from rdkit import Chem
from .graph import from_rdkit, fragment_package
from .inputs import from_document
from .provenance import check_scan_reference
from .preparation import prepare_geometry
from .gaussian import gaussian_input, molecular_charge, validate_level
from .api import get_backend, Q2Options, Q4Options, ScanOptions, resolve_options


def read_inputs(path):
    if path.suffix.lower() in {".sdf", ".sd"}:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        for i, rd in enumerate(supplier):
            if rd is None or not rd.GetNumConformers() or not rd.GetConformer().Is3D():
                raise ValueError(
                    f"SDF record {i}: invalid molecule or missing 3D coordinates"
                )
            if Chem.AddHs(rd).GetNumAtoms() != rd.GetNumAtoms():
                raise ValueError(f"SDF record {i}: explicit hydrogens required")
            name = rd.GetProp("_Name") if rd.HasProp("_Name") else str(i)
            m = from_rdkit(rd, name)
            for atom, pos in zip(m.Atoms, rd.GetConformer().GetPositions()):
                atom.coor = pos.tolist()
            types = (
                json.loads(rd.GetProp("QO2_atom_types"))
                if rd.HasProp("QO2_atom_types")
                else None
            )
            yield m, types
    else:
        with path.open() as stream:
            for line in stream:
                if not line.strip():
                    continue
                d = json.loads(line)
                if "coordinates" not in d:
                    raise ValueError("JSONL record requires coordinates in angstrom")
                yield from_document(d), d.get("atom_types")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["q0", "q1", "q2", "q4"])
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New directory; existing paths are refused",
    )
    parser.add_argument(
        "--level",
        required=True,
        help="Gaussian method/basis, e.g. B3LYP/def2SVP; no default",
    )
    parser.add_argument(
        "--multiplicity",
        type=int,
        required=True,
        help="Explicit spin multiplicity for every input record",
    )
    parser.add_argument(
        "--alkane-policy", choices=["supplied-types", "keep"], default="supplied-types"
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=12,
        help="Number of Q1 scan increments",
    )
    parser.add_argument(
        "--increment", type=float, default=30.0, help="Q1 increment in degrees"
    )
    minima_args = parser.add_mutually_exclusive_group()
    minima_args.add_argument(
        "--minima",
        type=Path,
        help="Q4 JSON: reference id -> torsion string -> minima angles (degrees)",
    )
    minima_args.add_argument(
        "--scan-records",
        type=Path,
        help="Q4 JSONL: legacy uses saved labels; qq detects minima from complete periodic scans",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=64,
        help="Q4 first N Cartesian-product combinations",
    )
    parser.add_argument(
        "--q4-optimization",
        choices=["none"],
        help="Required for Q4 CLI; explicit original optimize_structure=False branch",
    )
    parser.add_argument("--backend", choices=["legacy", "qq"], default="qq")
    parser.add_argument(
        "--bond-step", type=float, help="qq Q2 bond displacement in angstrom"
    )
    parser.add_argument(
        "--angle-step", type=float, help="qq Q2 angle displacement in degrees"
    )
    parser.add_argument(
        "--clash-scale",
        type=float,
        help="qq Q2/Q4 contact radius scale; zero disables contact rejection",
    )
    parser.add_argument("--seed", type=int, help="qq Q4 reproducible sampling seed")
    parser.add_argument(
        "--rmsd-threshold",
        type=float,
        help="qq Q4 duplicate RMSD threshold in angstrom",
    )
    parser.add_argument(
        "--rmsd-atoms", choices=["all", "heavy"], help="qq Q4 comparison atoms"
    )
    parser.add_argument(
        "--scan-step",
        type=float,
        help="qq scan-record grid spacing in degrees; default 30",
    )
    parser.add_argument(
        "--scan-energy-tolerance",
        type=float,
        help="qq plateau tolerance in input energy units; default 1e-8",
    )
    parser.add_argument(
        "--empty-minima",
        choices=["error", "reference"],
        help="qq Q4 empty-angle-list policy; default error",
    )
    parser.add_argument(
        "--scan-error-policy",
        choices=["strict", "collect"],
        help="qq scan failure isolation; default strict",
    )
    args = parser.parse_args(argv)
    if (
        args.scan_step is not None
        or args.scan_energy_tolerance is not None
        or args.scan_error_policy is not None
    ) and not (args.backend == "qq" and args.stage == "q4" and args.scan_records):
        parser.error("Scan options require qq Q4 --scan-records")
    engine = get_backend(args.backend)
    available = {
        "q2": {"bond_step", "angle_step", "clash_scale"},
        "q4": {"seed", "rmsd_threshold", "rmsd_atoms", "clash_scale", "empty_minima"},
    }
    specified = {
        key: getattr(args, key)
        for key in {
            "bond_step",
            "angle_step",
            "clash_scale",
            "seed",
            "rmsd_threshold",
            "rmsd_atoms",
            "empty_minima",
        }
        if getattr(args, key) is not None
    }
    if specified and (
        args.backend != "qq" or set(specified) - available.get(args.stage, set())
    ):
        parser.error(
            "qq options require the qq backend and their corresponding Q2 or Q4 stage"
        )
    backend_arguments = {}
    if args.backend == "qq" and args.stage in {"q2", "q4"}:
        try:
            cls = Q2Options if args.stage == "q2" else Q4Options
            backend_arguments = resolve_options(
                args.backend, args.stage, cls(**specified)
            )
        except (TypeError, ValueError) as exc:
            parser.error(str(exc))
    if args.stage == "q4" and (
        not (args.minima or args.scan_records) or args.q4_optimization is None
    ):
        parser.error(
            "Q4 requires --minima or --scan-records and --q4-optimization none"
        )
    if args.stage != "q4" and (
        args.minima or args.scan_records or args.q4_optimization
    ):
        parser.error(
            "Minima, scan records and Q4 optimization options apply only to Q4"
        )
    if args.limit < 1:
        parser.error("--limit must be positive")
    minima = {}
    if args.minima:
        minima = json.loads(args.minima.read_text())
    scan_report = None
    scan_result = None
    if args.scan_records:
        rows = [
            json.loads(line)
            for line in args.scan_records.read_text().splitlines()
            if line.strip()
        ]
        try:
            if args.backend == "qq":
                scan_result = engine.minima_from_scan_records(
                    rows,
                    options=ScanOptions(
                        error_policy=args.scan_error_policy or "strict",
                        expected_step=(
                            args.scan_step if args.scan_step is not None else 30.0
                        ),
                        energy_tolerance=(
                            args.scan_energy_tolerance
                            if args.scan_energy_tolerance is not None
                            else 1e-8
                        ),
                    ),
                    return_report=True,
                )
                minima, scan_report = scan_result["minima"], scan_result["report"]
            else:
                minima = engine.minima_from_scan_records(rows)
        except (ValueError, TypeError, KeyError) as exc:
            parser.error(str(exc))
    if (
        args.multiplicity < 1
        or args.steps < 1
        or not math.isfinite(args.increment)
        or args.increment == 0
    ):
        parser.error("Invalid multiplicity/scan settings")
    try:
        validate_level(args.level)
    except ValueError as exc:
        parser.error(str(exc))
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(
        stage=args.stage.upper(),
        backend=args.backend,
        backend_settings=(
            asdict(backend_arguments["options"]) if backend_arguments else None
        ),
        q4_executed=args.stage == "q4",
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        parameters={
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        index_convention="JSON zero-based; Gaussian one-based",
        environment=dict(
            python=platform.python_version(),
            rdkit=importlib.metadata.version("rdkit"),
            numpy=np.__version__,
        ),
        pair_policy="graph shortest distance 3 for Pair14; greater than 3 for Pair1n; rebuilt explicitly",
        runtime_sha256={
            str(p.relative_to(root.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
            for root in [Path(__file__).parent, Path(fragment_package.__file__).parent]
            for p in sorted(root.rglob("*.py"))
        },
        records=[],
    )
    if scan_report is not None:
        report["scan_analysis"] = scan_report
    if args.minima or args.scan_records:
        report["minima_input_sha256"] = hashlib.sha256(
            (args.minima or args.scan_records).read_bytes()
        ).hexdigest()
    failed = bool(scan_report and scan_report.get("failed_molecules"))
    try:
        for idx, (raw, types) in enumerate(read_inputs(args.input)):
            record = dict(index=idx, id=raw.name, input_metadata=raw.input_metadata)
            report["records"].append(record)
            try:
                m = prepare_geometry(
                    raw,
                    types,
                    args.alkane_policy if args.stage in {"q1", "q2"} else "keep",
                    backend=args.backend,
                    stage=args.stage,
                )
                charge = molecular_charge(m, args.multiplicity)
                record.update(charge=charge, multiplicity=args.multiplicity)
                ignore = args.alkane_policy == "supplied-types"
                if args.stage == "q0":
                    candidates = {"reference": np.asarray([a.coor for a in m.Atoms])}
                elif args.stage == "q4":
                    if args.backend == "qq":
                        record["scan_reference_check"] = check_scan_reference(
                            m,
                            scan_result,
                            level=args.level,
                            multiplicity=args.multiplicity,
                        )
                    candidates, q4_report = engine.generate_q4(
                        m, minima[str(raw.name)], limit=args.limit, **backend_arguments
                    )
                    record.update(q4_report)
                elif args.stage == "q1":
                    if args.backend == "qq":
                        q1_report = engine.generate_q1(m, ignore, return_report=True)
                        q1_report["provenance"].update(
                            level=args.level,
                            charge=charge,
                            multiplicity=args.multiplicity,
                        )
                        record.update(q1_report)
                        torsions = q1_report["torsions"]
                    else:
                        torsions = engine.generate_q1(m, ignore)
                    candidates = {
                        f"torsion-{i+1}": np.asarray([a.coor for a in m.Atoms])
                        for i in range(len(torsions))
                    }
                    record["torsions"] = torsions
                else:
                    candidates, q2_report = engine.generate_q2(
                        m, ignore, **backend_arguments
                    )
                    record.update(q2_report)
                record["files"] = []
                for i, (name, coords) in enumerate(candidates.items()):
                    filename = f"{idx:06d}_{args.stage}_{name}.gjf"
                    text = gaussian_input(
                        m,
                        coords,
                        charge,
                        args.multiplicity,
                        args.level,
                        args.stage.upper(),
                        name,
                        torsion=torsions[i] if args.stage == "q1" else None,
                        steps=args.steps,
                        increment=args.increment,
                        backend=args.backend,
                    )
                    (args.output / filename).write_text(text)
                    record["files"].append(
                        dict(
                            file=filename,
                            sha256=hashlib.sha256(text.encode()).hexdigest(),
                        )
                    )
                record["accepted_count"] = len(candidates)
                record["status"] = "ok"
            except (
                ValueError,
                FloatingPointError,
                IndexError,
                StopIteration,
                KeyError,
                TypeError,
                RuntimeError,
            ) as exc:
                failed = True
                record.update(status="error", error=str(exc))
    except (ValueError, KeyError, TypeError, OSError, IndexError, RuntimeError) as exc:
        failed = True
        report["input_error"] = str(exc)
    if not report["records"]:
        failed = True
        report.setdefault("input_error", "No input records")
    report["status"] = "errors" if failed else "ok"
    (args.output / "manifest.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            dict(
                status=report["status"],
                records=len(report["records"]),
                manifest=str(args.output / "manifest.json"),
            )
        )
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
