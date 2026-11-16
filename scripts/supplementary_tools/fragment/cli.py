"""Molecular fragmentation from local molecular files."""

import argparse
from pathlib import Path
import sys
from .api import DEFAULT_BACKEND


import gzip
import hashlib
import json


def dump(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, default=str)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        for i, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise ValueError("Each line must be a JSON object")
                record.setdefault("id", str(i))
                yield i, record, None
            except Exception as exc:
                yield i, None, f"{type(exc).__name__}: {exc}"


def generate(args):
    if args.backend != "f147" and args.config is not None:
        raise ValueError("--config requires f147; use --qq-config for QQ")
    if args.backend == "qq":
        from .qq.cli import generate as run
    elif args.backend == "f147":
        from .modern.batch import generate as run
    else:
        raise ValueError("Unknown backend. See README.md")
    return run(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog="Usage: README.md")
    parser.add_argument(
        "--version", action="version", version=f"1.3.0 / default: {DEFAULT_BACKEND}"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    g = commands.add_parser(
        "generate", help="Generate partitions and backend-specific outputs"
    )
    g.add_argument(
        "--input",
        type=Path,
        required=True,
        help="f147: .smi/.smiles, .sdf/.mol, .json, .jsonl (also .gz); qq: JSONL",
    )
    g.add_argument("--output", type=Path, required=True)
    g.add_argument(
        "--families", help="Comma-separated families; defaults depend on backend"
    )
    g.add_argument(
        "--backend",
        default=DEFAULT_BACKEND,
        help="f147 (default) or qq",
    )
    g.add_argument("--config", type=Path, help="FragmentConfig JSON for f147")
    g.add_argument(
        "--qq-config",
        type=Path,
        help="JSON export limits: min_heavy_atoms, max_heavy_atoms, max_candidates",
    )
    g.add_argument(
        "--include-parent",
        action="store_true",
        help="Also export the parent identity (f147)",
    )
    g.add_argument(
        "--allow-errors",
        action="store_true",
        help="Permit partial batch success with exit code 0",
    )
    args = parser.parse_args()
    try:
        return generate(args)
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
