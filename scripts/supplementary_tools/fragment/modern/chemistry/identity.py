"""Exact graph identity; Standard InChI remains auxiliary, not the merge key."""

import hashlib
from rdkit import Chem
from .input import check_stereo


def identify(mol):
    check_stereo(mol)
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(clean)
    Chem.AssignStereochemistry(clean, cleanIt=True, force=True)
    smiles = Chem.MolToSmiles(Chem.RemoveHs(clean), canonical=True, isomericSmiles=True)
    value = dict(
        identity_key=hashlib.sha256(("f_graph_v1\0" + smiles).encode()).hexdigest(),
        identity_policy="f_graph_v1_canonical_isomeric_smiles",
        smiles=smiles,
        formal_charge=sum(a.GetFormalCharge() for a in clean.GetAtoms()),
        radical_electrons=sum(a.GetNumRadicalElectrons() for a in clean.GetAtoms()),
        elements=sorted({a.GetSymbol() for a in clean.GetAtoms()}),
        inchi=None,
        inchi_key=None,
    )
    try:
        inchi = Chem.MolToInchi(clean)
        key = Chem.InchiToInchiKey(inchi) if inchi else None
        if not key:
            raise ValueError("InChI generation returned no identity")
        value.update(inchi=inchi, inchi_key=key)
    except Exception as exc:
        value["auxiliary_identity_error"] = f"{type(exc).__name__}: {exc}"
    return value
