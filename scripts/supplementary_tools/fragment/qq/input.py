"""Read explicit graphs or SMILES without using historical fragment fields."""

from rdkit import Chem
from ..common.graph import from_document

SUPPORTED = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}


def check_stereo(mol):
    """v1 identity supports ordinary tetrahedral and double-bond stereo only."""
    if mol.GetStereoGroups():
        raise ValueError("Enhanced stereo groups are outside QQ v1 identity support")
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
        raise ValueError("Unsupported stereo type for QQ v1 identity")


def read_parent(record, config):
    if "elements" in record:
        graph = from_document(record)
        if "isotopes" in record:
            if len(record["isotopes"]) != len(graph.Atoms):
                raise ValueError("Inconsistent isotope field length")
            for atom, isotope in zip(graph.Atoms, record["isotopes"]):
                if type(isotope) is not int or isotope < 0:
                    raise ValueError("Invalid isotope mass number")
                atom.isotope = isotope
        mol = graph.to_rdkit()
        status = "explicit_graph_no_stereo_inference"
        # Graph coordinates are not a reliable specification of stereo. A supplied
        # SMILES may provide it, but do not arbitrarily choose an ambiguous map.
        if record.get("smiles"):
            parsed = Chem.MolFromSmiles(record["smiles"])
            if parsed is None:
                raise ValueError(
                    "Supplied SMILES is unreadable; graph not silently substituted"
                )
            reference = Chem.AddHs(parsed)
            check_stereo(reference)
            Chem.Kekulize(reference, clearAromaticFlags=True)
            query = Chem.Mol(mol)
            Chem.Kekulize(query, clearAromaticFlags=True)
            matches = reference.GetSubstructMatches(query, uniquify=False, maxMatches=2)
            if reference.GetNumAtoms() != mol.GetNumAtoms() or not matches:
                raise ValueError("SMILES conflicts with explicit graph")
            if "isotopes" in record:
                graph_check, reference_check = Chem.Mol(mol), Chem.Mol(reference)
                Chem.RemoveStereochemistry(graph_check)
                Chem.RemoveStereochemistry(reference_check)
                if Chem.MolToSmiles(graph_check) != Chem.MolToSmiles(reference_check):
                    raise ValueError("SMILES conflicts with explicit isotope graph")
            sensitive = any(
                a.GetIsotope() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                for a in reference.GetAtoms()
            ) or any(
                b.GetStereo() != Chem.BondStereo.STEREONONE
                for b in reference.GetBonds()
            )
            # Do not choose an arbitrary mapping or silently discard identity-
            # bearing information. Users can supply SMILES or an explicit Mol.
            if len(matches) > 1 and sensitive:
                raise ValueError(
                    "Ambiguous SMILES-to-graph identity mapping; supply SMILES alone or an explicit RDKit Mol"
                )
            if len(matches) == 1:
                mol = Chem.RenumberAtoms(reference, list(matches[0]))
                status = "stereo_from_unique_smiles_graph_mapping"
            else:
                status = "ambiguous_smiles_mapping_graph_only"
    elif "mol" in record:
        if not isinstance(record["mol"], Chem.Mol):
            raise ValueError("mol must be an RDKit molecule")
        mol = Chem.Mol(record["mol"])
        status = "explicit_rdkit_input"
    else:
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is None:
            raise ValueError("SMILES could not be parsed")
        status = "input_smiles"
    if not mol.GetNumAtoms():
        raise ValueError("Empty molecule")
    Chem.SanitizeMol(mol)
    check_stereo(mol)
    unknown = sorted(
        {a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() not in SUPPORTED}
    )
    if unknown:
        raise ValueError("Unsupported elements: " + ", ".join(unknown))
    if not config.allow_radicals and any(
        a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        raise ValueError("Radical input requires allow_radicals=True")
    mol = Chem.AddHs(mol)
    mol.RemoveAllConformers()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol, status
