"""Validate requests and schedule parent processing; persistence lives in io."""

from collections import Counter
import json
import time
from .config import FragmentConfig
from .pipeline import process, validate_request
from .io.jsonl import digest
from .io.inputs import input_format, records
from .io.store import BatchStore
from .io.manifest import write_manifest


def generate(args):
    if args.qq_config is not None:
        raise ValueError("The f147 backend does not accept --qq-config")
    values = json.loads(args.config.read_text()) if args.config else {}
    if not isinstance(values, dict):
        raise ValueError("--config must contain a JSON object")
    config = FragmentConfig(**values)
    families = None if args.families is None else tuple(args.families.split(","))
    selected = validate_request(args.backend, families, args.include_parent, config)
    input_format(args.input)
    input_hash = digest(args.input)
    counts = Counter()
    started = time.monotonic()
    with BatchStore(args.output) as store:
        for line, record, error in records(args.input):
            counts["input_records"] += 1
            if error is not None:
                counts["input_failures"] += 1
                store.write_error(dict(input_line=line, stage="input", error=error))
                continue
            try:
                result = process(
                    record,
                    families,
                    args.include_parent,
                    backend=args.backend,
                    config=config,
                )
            except Exception as exc:
                counts["input_failures"] += 1
                store.write_error(
                    dict(
                        input_line=line,
                        input_id=record.get("id"),
                        stage="parent",
                        error=f"{type(exc).__name__}: {exc}",
                    )
                )
                continue
            counts["processed_parents"] += 1
            counts["elementary_blocks"] += len(result["elementary_fragments"])
            # Preserve sparse counts in the existing output schema.
            for key, field in (
                ("stage_errors", "errors"),
                ("excluded", "excluded"),
                ("occurrences", "occurrences"),
            ):
                if result[field]:
                    counts[key] += len(result[field])
            store.write_parent(result, line)
            if counts["processed_parents"] % 100 == 0:
                store.commit()
        unique = store.export()
        if unique:
            counts["unique_identities"] += unique
    write_manifest(
        args.output,
        backend=args.backend,
        selected=selected,
        include_parent=args.include_parent,
        input_hash=input_hash,
        config=config,
        counts=counts,
        elapsed_seconds=time.monotonic() - started,
    )
    print(json.dumps(dict(counts), indent=2))
    return (
        2
        if (counts["input_failures"] or counts["stage_errors"])
        and not args.allow_errors
        else 0
    )
