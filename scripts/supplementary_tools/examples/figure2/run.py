"""Run one modern fragmentation example and draw its coloured Figure 2."""

import argparse
from pathlib import Path
import sys
import shutil
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
sys.path.insert(0, str(PACKAGE))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=HERE / "input" / "figure2_3702543.smi",
        help="One parent in any modern input format (default: figure2_3702543).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "output",
        help="New result directory (default: beside this script).",
    )
    parser.add_argument("--backend", choices=("f147", "f390"), default="f147")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output directory after explicit confirmation.",
    )
    args = parser.parse_args()
    # Check drawing dependencies before creating any output.
    try:
        from drawing import draw_figure
        from fragment.modern.io.inputs import records
        from fragment.modern.batch import generate
        import json

        iterator = records(args.input)
        try:
            first = next(iterator, None)
            if first is None or next(iterator, None) is not None:
                raise ValueError(
                    "This figure example requires exactly one parent record."
                )
        finally:
            iterator.close()
        position, record, error = first
        if error:
            raise ValueError(f"Input record {position}: {error}")
        if args.output.exists():
            if not args.overwrite:
                raise ValueError(
                    f"Output already exists: {args.output}; use --overwrite to replace it."
                )
            if not args.output.is_dir():
                raise ValueError(f"Output path is not a directory: {args.output}")
            shutil.rmtree(args.output)
        status = generate(
            SimpleNamespace(
                input=args.input,
                output=args.output,
                backend=args.backend,
                config=None,
                qq_config=None,
                descriptors_only=False,
                families=None,
                include_parent=False,
                allow_errors=False,
            )
        )
        if status:
            print(
                f'Fragmentation failed; inspect {args.output / "errors.jsonl"}. No figure created.',
                file=sys.stderr,
            )
            return status
        parent = json.loads((args.output / "parents.jsonl").read_text())
        occurrences = [
            json.loads(line)
            for line in (args.output / "occurrences.jsonl").read_text().splitlines()
        ]
        draw_figure(record, parent, occurrences, args.output)
        print(
            f'\nResults: {args.output.resolve()}\nFigure: {(args.output / "figure2.pdf").resolve()}'
        )
        return 0
    except ImportError as exc:
        print(
            f'{exc}\nInstall package dependencies: python -m pip install -r {PACKAGE / "requirements.txt"}',
            file=sys.stderr,
        )
        return 1
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
