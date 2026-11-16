"""147 partitioning with per-call state."""

from rdkit import Chem

from .grouping import Grouping
from .connectivity import find_mole
from ..chemistry.preprocess import prepare


class _Bond:

    def __init__(self, a, b):
        self.a1, self.a2 = (a, b)

    def get_type(self, molecule):
        atom = molecule.Atoms[self.a1]
        return atom.bond_type[atom.connect.index(self.a2)]


class _RuntimeMolScalpel(Grouping):
    """Guard residual cycles without altering the archived grouping kernel."""

    def run(self, molecule):
        # Structural ranks resolve equal-priority chain directions and terminal
        # merges without using input atom numbers. Preserve all parent indices.
        ranked = Chem.Mol(molecule.rdmol)
        Chem.SetAromaticity(ranked)
        molecule._stable_atom_ranks = tuple(
            Chem.CanonicalRankAtoms(ranked, breakTies=False, includeAtomMaps=False)
        )
        return super().run(molecule)

    def _molecule_split(self, molecule):
        # Late heavy-atom merges can leave attached H with an obsolete group.
        # Synchronize ownership before computing boundary flags and descriptors.
        # Isotopes stay on their original atoms; no hydrogens are regenerated.
        for atom in molecule.Atoms:
            if atom.element != "H" or len(atom.connect) != 1:
                continue
            neighbor = molecule.Atoms[atom.connect[0]]
            if neighbor.element == "H" or atom.bond_type[0] != "1":
                continue
            if atom.function_group_id != neighbor.function_group_id:
                atom.function_group = neighbor.function_group
                atom.function_group_id = neighbor.function_group_id
                atom.function_group_tag = "H_alpha"
        return super()._molecule_split(molecule)

    def _ane_chain(self, molecule, res_atom):
        members = set(res_atom)
        adjacency = {
            a: [b for b in molecule.Atoms[a].connect if b in members] for a in res_atom
        }
        remaining = []
        for component in find_mole(adjacency):
            edges = sum((len(adjacency[a]) for a in component)) // 2
            if edges < len(component):
                remaining.extend(component)
                continue
            group_ids = {molecule.Atoms[a].function_group_id for a in component}
            new_id = min(group_ids)
            for atom in molecule.Atoms:
                if atom.function_group_id in group_ids:
                    atom.function_group_id = new_id
            molecule.runtime_adjustments = getattr(molecule, "runtime_adjustments", [])
            molecule.runtime_adjustments.append(
                {"kind": "retained_residual_cycle", "atoms": sorted(component)}
            )
        return super()._ane_chain(molecule, remaining) if remaining else []


def partition(molecule):
    """Run historical typing and partitioning, retaining boundary flags."""
    supported = {
        "C",
        "H",
        "N",
        "O",
        "S",
        "P",
        "F",
        "Cl",
        "Br",
        "I",
        "Si",
        "B",
        "Na",
        "K",
        "Ca",
        "Mg",
    }
    unsupported = sorted(set(molecule.elements) - supported)
    if unsupported:
        raise ValueError(
            "f147 has no typing rules for element(s): " + ", ".join(unsupported)
        )
    molecule = prepare(molecule)
    molecule._rings = molecule.rings
    molecule.Bonds = [_Bond(b.a1, b.a2) for b in molecule.Bonds]
    _RuntimeMolScalpel().run(molecule)
    return molecule
