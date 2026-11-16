"""Explicit, immutable settings for the new qq algorithms (not DFT-calibrated)."""

from dataclasses import dataclass
import math


def _finite(name, value, lower, upper=None):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ValueError(f"{name} must be a finite number")
    if value < lower or (upper is not None and value > upper):
        raise ValueError(f"{name} is outside its allowed range")


@dataclass(frozen=True)
class Q2Options:
    bond_step: float = 0.05
    angle_step: float = 5.0
    clash_scale: float = 0.6

    def __post_init__(self):
        _finite("bond_step", self.bond_step, 1e-8)
        _finite("angle_step", self.angle_step, 1e-8, 179.0)
        _finite("clash_scale", self.clash_scale, 0.0, 1.0)


@dataclass(frozen=True)
class Q4Options:
    seed: int = 0
    rmsd_threshold: float = 0.1
    rmsd_atoms: str = "all"
    clash_scale: float = 0.6
    empty_minima: str = "error"

    def __post_init__(self):
        if self.empty_minima not in {"error", "reference"}:
            raise ValueError("empty_minima must be error or reference")
        if type(self.seed) is not int:
            raise ValueError("seed must be an integer")
        _finite("rmsd_threshold", self.rmsd_threshold, 0.0)
        _finite("clash_scale", self.clash_scale, 0.0, 1.0)
        if self.rmsd_atoms not in {"all", "heavy"}:
            raise ValueError("rmsd_atoms must be all or heavy")


@dataclass(frozen=True)
class ScanOptions:
    expected_step: float = 30.0
    energy_tolerance: float = 1e-8
    angle_tolerance: float = 1e-6
    error_policy: str = "strict"

    def __post_init__(self):
        if self.error_policy not in {"strict", "collect"}:
            raise ValueError("error_policy must be strict or collect")
        _finite("expected_step", self.expected_step, 1e-6, 120.0)
        _finite("energy_tolerance", self.energy_tolerance, 0.0)
        _finite("angle_tolerance", self.angle_tolerance, 1e-12)
        if self.angle_tolerance >= self.expected_step / 2:
            raise ValueError("angle_tolerance must be less than half the scan step")
        if abs(360 / self.expected_step - round(360 / self.expected_step)) > 1e-8:
            raise ValueError("expected_step must divide 360 degrees")
