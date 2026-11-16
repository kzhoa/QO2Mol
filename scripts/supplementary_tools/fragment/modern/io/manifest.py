"""Write the versioned run manifest and runtime fingerprint."""

from pathlib import Path
import hashlib
import platform
from rdkit import rdBase
from .jsonl import dump


def write_manifest(
    output,
    *,
    backend,
    selected,
    include_parent,
    input_hash,
    config,
    counts,
    elapsed_seconds
):
    root = Path(__file__).resolve().parents[1]
    manifest = dict(
        backend=backend,
        profile=backend + "_v1",
        tool_version="1.2.0",
        families=list(selected),
        include_parent=include_parent,
        input_sha256=input_hash,
        config=config.to_dict(),
        config_sha256=config.fingerprint,
        counts=dict(counts),
        python=platform.python_version(),
        rdkit=rdBase.rdkitVersion,
        selection="all_successful_exact_graph_identities",
        historical_reproduction=False,
        product_normalization="rdkit_memory_boundary_v1",
        coordinates_generated=False,
        runtime_sha256={
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix in (".py", ".json", ".def", ".yaml")
        },
        status="errors" if counts["input_failures"] or counts["stage_errors"] else "ok",
        elapsed_seconds=elapsed_seconds,
    )
    (output / "run.json").write_text(dump(manifest) + "\n")
