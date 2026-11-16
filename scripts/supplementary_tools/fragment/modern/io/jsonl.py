"""JSON Lines and input digest helpers."""

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
    opener = gzip.open if path.suffix.lower() == ".gz" else open
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
