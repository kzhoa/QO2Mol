"""Validated resource/export policy, independent of the fragmentation kernel."""

from dataclasses import asdict, dataclass
import hashlib
import json


@dataclass(frozen=True)
class FragmentConfig:
    min_heavy_atoms: int = 1
    max_heavy_atoms: int | None = None
    max_input_atoms: int = 10000
    max_products_per_parent: int = 10000
    allow_radicals: bool = False

    def __post_init__(self):
        for name in ("min_heavy_atoms", "max_input_atoms", "max_products_per_parent"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(name + " must be a positive integer")
        if self.max_heavy_atoms is not None:
            if (
                type(self.max_heavy_atoms) is not int
                or self.max_heavy_atoms < self.min_heavy_atoms
            ):
                raise ValueError(
                    "max_heavy_atoms must be an integer >= min_heavy_atoms, or null"
                )
        if type(self.allow_radicals) is not bool:
            raise TypeError("allow_radicals must be boolean")

    def to_dict(self):
        return asdict(self)

    @property
    def fingerprint(self):
        return hashlib.sha256(
            json.dumps(self.to_dict(), sort_keys=True).encode()
        ).hexdigest()
