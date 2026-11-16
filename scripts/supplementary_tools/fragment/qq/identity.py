"""Exact graph identity independent of the auxiliary Standard InChI identity."""

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
    result = {
        "identity_key": hashlib.sha256(
            ("qq_strict_v1\0" + smiles).encode()
        ).hexdigest(),
        "identity_policy": "qq_strict_v1_canonical_isomeric_smiles",
        "smiles": smiles,
        "formal_charge": sum(a.GetFormalCharge() for a in clean.GetAtoms()),
        "radical_electrons": sum(a.GetNumRadicalElectrons() for a in clean.GetAtoms()),
        "elements": sorted({a.GetSymbol() for a in clean.GetAtoms()}),
        "inchi": None,
        "inchi_key": None,
    }
    try:
        inchi = Chem.MolToInchi(clean)
        if not inchi:
            raise ValueError("InChI generation returned no identity")
        result.update(inchi=inchi, inchi_key=Chem.InchiToInchiKey(inchi))
    except Exception as exc:
        result["auxiliary_identity_error"] = f"{type(exc).__name__}: {exc}"
    return result
