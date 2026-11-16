"""Read modern molecular inputs into the shared process() record contract."""

import gzip
import json
from pathlib import Path
from rdkit import Chem
from .jsonl import records as jsonl_records

FORMATS = {
    ".jsonl": "jsonl",
    ".ndjson": "jsonl",
    ".json": "json",
    ".smi": "smiles",
    ".smiles": "smiles",
    ".sdf": "sdf",
    ".mol": "mol",
}


def input_format(path):
    path = Path(path)
    suffix = Path(path.stem).suffix if path.suffix.lower() == ".gz" else path.suffix
    try:
        return FORMATS[suffix.lower()]
    except KeyError:
        raise ValueError(
            "Unsupported input extension; use .smi/.smiles, .sdf/.mol, .json or .jsonl (optionally .gz)"
        ) from None


def _object(value, index):
    if not isinstance(value, dict):
        raise ValueError("Each record must be a JSON object")
    value.setdefault("id", str(index))
    return value


def _molecule(mol, index):
    if mol is None:
        raise ValueError("Could not parse molecular structure")
    name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    record = {"id": name or str(index), "mol": mol}
    # Only copy metadata, never let SD properties override molecular structure.
    for key in ("id", "source", "collection", "scope"):
        if mol.HasProp(key) and mol.GetProp(key).strip():
            record[key] = mol.GetProp(key)
    return record


def records(path):
    """Yield (position, record, error); positions are 1-based lines or records.

    JSONL/SMILES use physical line numbers; JSON arrays and SDF use record
    ordinals. Single JSON objects and MOL files have position 1.
    """
    path = Path(path)
    kind = input_format(path)
    if kind == "jsonl":
        yield from jsonl_records(path)
        return
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    if kind in ("json", "smiles"):
        with opener(path, "rt", encoding="utf-8-sig") as stream:
            if kind == "json":
                try:
                    data = json.load(stream)
                except (ValueError, UnicodeError) as exc:
                    yield 1, None, f"{type(exc).__name__}: {exc}"
                    return
                rows = enumerate(data if isinstance(data, list) else [data], 1)
            else:
                rows = enumerate(stream, 1)
            for index, value in rows:
                try:
                    if kind == "smiles":
                        value = value.strip()
                        if not value or value.startswith("#"):
                            continue
                        fields = value.split(maxsplit=1)
                        record = {
                            "smiles": fields[0],
                            "id": fields[1] if len(fields) > 1 else str(index),
                        }
                    else:
                        record = _object(value, index)
                except ValueError as exc:
                    yield index, None, f"{type(exc).__name__}: {exc}"
                else:
                    yield index, record, None
        return
    with opener(path, "rb") as stream:
        if kind == "sdf":
            molecules = Chem.ForwardSDMolSupplier(
                stream, sanitize=True, removeHs=False, strictParsing=True
            )
        else:
            try:
                mol = Chem.MolFromMolBlock(
                    stream.read().decode("utf-8-sig"),
                    sanitize=True,
                    removeHs=False,
                    strictParsing=True,
                )
            except (ValueError, RuntimeError, UnicodeError) as exc:
                yield 1, None, f"{type(exc).__name__}: {exc}"
                return
            molecules = [mol]
        for index, mol in enumerate(molecules, 1):
            try:
                record = _molecule(mol, index)
            except ValueError as exc:
                yield index, None, f"{type(exc).__name__}: {exc}"
            else:
                yield index, record, None
