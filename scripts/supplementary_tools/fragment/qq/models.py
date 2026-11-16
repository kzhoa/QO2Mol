"""Versioned configuration and graph descriptors for QQ bond-based fragmentation."""

from dataclasses import asdict, dataclass
import hashlib
import json
from functools import lru_cache
from pathlib import Path

VERSION = "qq_basic_units_v2"
TOOL_VERSION = "1.3.0"
FAMILIES = ("elementary",)


@dataclass(frozen=True)
class QQConfig:
    """Export/resource limits only; chemistry is the frozen basic-units rule set.

    Removed options (context, linker, sampling and rule overrides) are rejected
    by the constructor, including old options set to their neutral values.
    """
    min_heavy_atoms: int = 1
    max_heavy_atoms: int = 60
    max_candidates: int = 10000

    def __post_init__(self):
        for key in ("min_heavy_atoms", "max_heavy_atoms", "max_candidates"):
            if type(getattr(self, key)) is not int or getattr(self, key) < 1:
                raise ValueError(f"{key} must be a positive integer")
        if self.min_heavy_atoms > self.max_heavy_atoms:
            raise ValueError("min_heavy_atoms exceeds max_heavy_atoms")

    # Private rule stages share these read-only policy values with the frozen
    # implementation. They are not constructor options or configurable knobs.
    protect_conjugated = property(lambda self: True)
    protect_charged = property(lambda self: True)
    protect_terminal = property(lambda self: True)
    cut_acyclic_cc = property(lambda self: False)
    allow_radicals = property(lambda self: False)
    protected_smarts = property(lambda self: ())
    ether_linker_max_heavy = property(lambda self: 0)

    def to_dict(self):
        return asdict(self)

    @property
    def fingerprint(self):
        text = json.dumps({"version": VERSION, "config": self.to_dict()}, sort_keys=True)
        return hashlib.sha256(text.encode()).hexdigest()


@dataclass(frozen=True)
class BondDecision:
    bond_index: int
    atoms: tuple[int, int]
    cut: bool
    reason: str
    protections: tuple[str, ...]
    matched_cut_rules: tuple[str, ...]


@dataclass(frozen=True)
class ElementaryFragment:
    index: int
    atoms: tuple[int, ...]
    heavy_atoms: tuple[int, ...]
    contains_ring: bool


@lru_cache(maxsize=1)
def implementation_hash():
    """Fingerprint this core and its shared graph adapter, without local paths."""
    root = Path(__file__).resolve().parent
    files = sorted(root.glob("*.py")) + [root.parent / "common" / "graph.py"]
    digest = hashlib.sha256()
    for path in files:
        digest.update(path.name.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()
