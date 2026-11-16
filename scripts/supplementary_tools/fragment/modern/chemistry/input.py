"""Explicit molecular input with no coordinate-derived identity guesses."""

from rdkit import Chem
from .graph import from_document


def check_stereo(mol):
    if mol.GetStereoGroups():
        raise ValueError(
            "Enhanced stereo groups are not supported; supply explicit absolute stereo"
        )
    atom_types = {
        Chem.ChiralType.CHI_UNSPECIFIED,
        Chem.ChiralType.CHI_TETRAHEDRAL_CW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    }
    bond_types = {
        Chem.BondStereo.STEREONONE,
        Chem.BondStereo.STEREOANY,
        Chem.BondStereo.STEREOE,
        Chem.BondStereo.STEREOZ,
        Chem.BondStereo.STEREOCIS,
        Chem.BondStereo.STEREOTRANS,
    }
    if any(a.GetChiralTag() not in atom_types for a in mol.GetAtoms()) or any(
        b.GetStereo() not in bond_types for b in mol.GetBonds()
    ):
        raise ValueError("Unsupported stereochemical representation")


def read_parent(record, config):
    if not isinstance(record, dict):
        raise TypeError("A parent record must be a dictionary")
    if "elements" in record:
        graph = from_document(record)
        isotopes = record.get("isotopes", [0] * len(graph.Atoms))
        if len(isotopes) != len(graph.Atoms) or any(
            type(i) is not int or i < 0 for i in isotopes
        ):
            raise ValueError("Invalid isotope array")
        for atom, isotope in zip(graph.Atoms, isotopes):
            atom.isotope = isotope
        mol = graph.to_rdkit()
        status = "explicit_graph_no_stereo_inference"
        if record.get("smiles"):
            reference = Chem.MolFromSmiles(record["smiles"])
            if reference is None:
                raise ValueError("Invalid SMILES accompanying explicit graph")
            reference = Chem.AddHs(reference)
            check_stereo(reference)
            Chem.Kekulize(reference, clearAromaticFlags=True)
            query = Chem.Mol(mol)
            Chem.Kekulize(query, clearAromaticFlags=True)
            matches = reference.GetSubstructMatches(query, uniquify=False, maxMatches=2)
            if reference.GetNumAtoms() != mol.GetNumAtoms() or not matches:
                raise ValueError("SMILES conflicts with explicit graph")
            if "isotopes" in record:
                left, right = Chem.Mol(mol), Chem.Mol(reference)
                Chem.RemoveStereochemistry(left)
                Chem.RemoveStereochemistry(right)
                if Chem.MolToSmiles(left) != Chem.MolToSmiles(right):
                    raise ValueError("SMILES conflicts with explicit isotope graph")
            sensitive = any(
                a.GetIsotope() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                for a in reference.GetAtoms()
            ) or any(
                b.GetStereo() != Chem.BondStereo.STEREONONE
                for b in reference.GetBonds()
            )
            if len(matches) > 1 and sensitive:
                raise ValueError(
                    "Ambiguous identity mapping; supply SMILES alone or an indexed RDKit Mol"
                )
            if len(matches) == 1:
                mol = Chem.RenumberAtoms(reference, list(matches[0]))
                status = "unique_smiles_graph_mapping"
    elif "mol" in record:
        if not isinstance(record["mol"], Chem.Mol):
            raise TypeError("mol must be an RDKit Mol")
        mol, status = Chem.Mol(record["mol"]), "rdkit_molecule"
    else:
        if not isinstance(record.get("smiles"), str):
            raise ValueError("Provide smiles, mol, or an explicit molecular graph")
        mol, status = Chem.MolFromSmiles(record["smiles"]), "input_smiles"
        if mol is None:
            raise ValueError("SMILES could not be parsed")
    if not mol.GetNumAtoms():
        raise ValueError("Empty molecule")
    Chem.SanitizeMol(mol)
    check_stereo(mol)
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Dummy atoms are not supported")
    if not config.allow_radicals and any(
        a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        raise ValueError("Radical input requires allow_radicals=True")
    if any(
        b.GetBondType()
        not in (
            Chem.BondType.SINGLE,
            Chem.BondType.DOUBLE,
            Chem.BondType.TRIPLE,
            Chem.BondType.AROMATIC,
        )
        for b in mol.GetBonds()
    ):
        raise ValueError(
            "Only covalent single, double, triple and aromatic bonds are supported"
        )
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() > config.max_input_atoms:
        raise ValueError("max_input_atoms exceeded")
    mol.RemoveAllConformers()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol, status
