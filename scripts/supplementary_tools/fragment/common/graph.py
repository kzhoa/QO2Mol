"""Minimal graph adapter. No coordinates are interpreted as stereochemistry."""

from itertools import combinations
from types import SimpleNamespace
from rdkit import Chem


class Atom:
    def __init__(self, index):
        self.ID = index
        self.elem = "H"
        self.connect = []
        self.bond_type = []
        self.formal_charge = 0
        self.primitive_formal_charge = 0
        self.coor = [0.0, 0.0, 0.0]
        self.has_ring = []
        self.has_ring_size = []
        self.has_ring_property = []
        self.connectivity_type = []
        self.residu = ""
        self.atom_type_name = ""
        self.isotope = 0

    No = property(lambda s: s.ID)
    element = property(lambda s: s.elem, lambda s, v: setattr(s, "elem", v))
    symbol = property(lambda s: s.elem)
    connectivity = property(lambda s: s.connect, lambda s, v: setattr(s, "connect", v))
    inring = property(lambda s: s.has_ring, lambda s, v: setattr(s, "has_ring", v))
    ring_size = property(
        lambda s: s.has_ring_size, lambda s, v: setattr(s, "has_ring_size", v)
    )
    ring_prop = property(
        lambda s: s.has_ring_property, lambda s, v: setattr(s, "has_ring_property", v)
    )
    # Alias bond-type details to the corresponding connectivity-type field.
    bond_type_detail = property(
        lambda s: s.connectivity_type, lambda s, v: setattr(s, "connectivity_type", v)
    )


class Molecule:
    def __init__(self, name="molecule"):
        self.name = name
        self.Atoms = []
        self.steps = []
        self.rings = {}
        self.ring_blocks = []
        self.rdmol = None
        self.input_metadata = {}

    def create_atoms(self, n):
        self.Atoms = [Atom(i) for i in range(n)]

    def add_atom(self):
        self.Atoms.append(Atom(len(self.Atoms)))

    @property
    def elements(self):
        return [a.elem for a in self.Atoms]

    @property
    def connectivity(self):
        return [a.connect for a in self.Atoms]

    ring_dict = property(lambda s: s.rings)
    _ring_dict = property(lambda s: s.rings)
    ring_stru = property(lambda s: s.ring_blocks)

    def create_topols(self):
        self.Bonds = [
            SimpleNamespace(a1=a.ID, a2=j)
            for a in self.Atoms
            for j in a.connect
            if a.ID < j
        ]
        self.Angles = [
            SimpleNamespace(a1=i, a2=a.ID, a3=j)
            for a in self.Atoms
            for i, j in combinations(a.connect, 2)
        ]

    def to_rdkit(self):
        rw = Chem.RWMol()
        for a in self.Atoms:
            atom = Chem.Atom(a.elem)
            atom.SetFormalCharge(int(a.formal_charge))
            atom.SetIsotope(a.isotope)
            atom.SetNoImplicit(True)
            rw.AddAtom(atom)
        bt = {
            "1": Chem.BondType.SINGLE,
            "2": Chem.BondType.DOUBLE,
            "3": Chem.BondType.TRIPLE,
            "ar": Chem.BondType.AROMATIC,
        }
        for a in self.Atoms:
            for j, b in zip(a.connect, a.bond_type):
                if j > a.ID:
                    rw.AddBond(a.ID, j, bt[b])
        m = rw.GetMol()
        Chem.SanitizeMol(m)
        return m


def from_rdkit(rd, name="molecule"):
    rd = Chem.AddHs(Chem.Mol(rd))
    Chem.Kekulize(rd, clearAromaticFlags=True)
    m = Molecule(name)
    m.create_atoms(rd.GetNumAtoms())
    m.rdmol = rd
    for atom in rd.GetAtoms():
        a = m.Atoms[atom.GetIdx()]
        a.elem = atom.GetSymbol()
        a.formal_charge = atom.GetFormalCharge()
        a.primitive_formal_charge = a.formal_charge
        a.isotope = atom.GetIsotope()
        a.connect = [x.GetIdx() for x in atom.GetNeighbors()]
        a.bond_type = [
            str(int(rd.GetBondBetweenAtoms(a.ID, j).GetBondTypeAsDouble()))
            for j in a.connect
        ]
    m.create_topols()
    return m


def from_smiles(smiles, name="molecule"):
    rd = Chem.MolFromSmiles(smiles)
    if rd is None:
        raise ValueError("SMILES could not be parsed")
    return from_rdkit(rd, name)


def from_document(d):
    """Keep the stored atom indices and connectivity order exactly."""
    m = Molecule(str(d.get("id", d.get("inchi_key", "molecule"))))
    m.create_atoms(len(d["elements"]))
    n = len(m.Atoms)
    for key in ["connectivity", "bond_type", "formal_charge"]:
        if len(d[key]) != n:
            raise ValueError("Inconsistent graph field lengths: " + key)
    for i, a in enumerate(m.Atoms):
        a.elem = d["elements"][i]
        a.connect = [int(j) for j in d["connectivity"][i]]
        a.bond_type = [str(b) for b in d["bond_type"][i]]
        a.formal_charge = int(d["formal_charge"][i])
        a.primitive_formal_charge = a.formal_charge
        a.coor = list(d.get("coordinates", [[0.0, 0.0, 0.0]] * n)[i])
        if len(a.connect) != len(a.bond_type) or len(set(a.connect)) != len(a.connect):
            raise ValueError("Invalid connectivity row")
        for j, b in zip(a.connect, a.bond_type):
            if j < 0 or j >= n or j == i or i not in d["connectivity"][j]:
                raise ValueError("Asymmetric or invalid bond")
            if str(d["bond_type"][j][list(d["connectivity"][j]).index(i)]) != b:
                raise ValueError("Asymmetric bond order")
    m.create_topols()
    m.input_metadata = {
        k: d[k]
        for k in ["inchi_key", "inchi", "smiles", "frag_type", "source", "collection"]
        if k in d
    }
    return m


def attach_identity_graph(m):
    """Map stored SMILES onto the explicit graph, retaining available stereo/isotopes.

    A connectivity-only document cannot establish missing stereochemistry. The
    status is exported; cached historical identifiers are never overwritten.
    """
    rd = m.to_rdkit()
    Chem.Kekulize(rd, clearAromaticFlags=True)
    meta = m.input_metadata
    status = "graph_only_no_stereo"
    if meta.get("smiles"):
        parent = from_smiles(meta["smiles"]).rdmol
        query = Chem.Mol(rd)
        Chem.SanitizeMol(parent)
        Chem.SanitizeMol(query)
        # Saved structural fields lack isotope annotations.
        for a in query.GetAtoms():
            a.SetIsotope(0)
        matches = parent.GetSubstructMatches(query, uniquify=False, maxMatches=2)
        if parent.GetNumAtoms() != rd.GetNumAtoms() or not matches:
            raise ValueError("Stored SMILES does not match explicit graph")
        rd = Chem.RenumberAtoms(parent, list(matches[0]))
        status = "stereo_from_smiles_graph_mapping"
        if len(matches) > 1:
            status += ";symmetry_mapping_ambiguous"
    m.rdmol = rd
    return status
