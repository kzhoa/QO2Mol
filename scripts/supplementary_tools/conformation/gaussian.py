"""Shared Gaussian input formatting and electronic-state parameter checks."""

from .graph import Element


def validate_level(level):
    """Require one method/basis token, consistently for Python and CLI callers."""
    if (
        not isinstance(level, str)
        or any(c.isspace() for c in level)
        or "/" not in level
        or level.startswith(("#", "%"))
    ):
        raise ValueError("Expected a single method/basis token")


def molecular_charge(molecule, multiplicity):
    """Return formal net charge after checking spin/electron-count compatibility."""
    charge = sum(a.formal_charge for a in molecule.Atoms)
    electrons = sum(Element.get(a.elem).number for a in molecule.Atoms) - charge
    if (
        type(multiplicity) is not int
        or multiplicity < 1
        or electrons < multiplicity - 1
        or (electrons - multiplicity + 1) % 2
    ):
        raise ValueError("Multiplicity inconsistent with electron count")
    return charge


def gaussian_input(
    m,
    coordinates,
    charge,
    multiplicity,
    level,
    stage,
    name,
    torsion=None,
    steps=12,
    increment=30.0,
    backend="qq",
):
    job = {
        "Q0": "opt freq pop=chelpg",
        "Q1": "opt=modredundant",
        "Q2": "force",
        "Q4": "opt freq pop=chelpg",
    }[stage]
    if backend not in {"legacy", "qq"}:
        raise ValueError("backend must be legacy or qq")
    label = backend
    text = f"#p {level} {job} nosymm\n\nQO2Mol {label} {stage} {name}\n\n{charge} {multiplicity}\n"
    for atom, pos in zip(m.Atoms, coordinates):
        text += f"{atom.elem:2s} " + " ".join(f"{v:.10f}" for v in pos) + "\n"
    text += "\n"
    if torsion is not None:
        text += (
            "D "
            + " ".join(str(i + 1) for i in torsion)
            + f" S {steps} {increment:g}\n\n"
        )
    return text
