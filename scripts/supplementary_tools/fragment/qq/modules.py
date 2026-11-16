"""Attachment-aware assembly identities, separate from capped molecule identity."""

import hashlib
from collections import Counter
from rdkit import Chem

MODULE_POLICY = "qq_attachment_v1"


def describe_module(product, construction):
    """Describe the ports of a product constructed by :func:`products.cap`.

    Canonical atom maps encode 2**number_of_single_bond_ports on each retained
    heavy atom (1 means no port). Original atom maps and parent indices never
    enter the identity. The capped graph supplies isotope, charge and surviving
    stereo information. Cut stereochemistry is not reconstructed by this key.
    This is an assembly descriptor, not a different molecule for DFT.
    """
    mapping = construction["product_to_parent"]
    if len(mapping) != product.GetNumAtoms():
        raise ValueError("Product mapping length mismatch")
    inverse = {old: new for new, old in enumerate(mapping) if old is not None}
    ports = []
    counts = Counter()
    for edge in construction["boundary_bonds"]:
        old = edge["kept_atom"]
        new = inverse[old]
        if product.GetAtomWithIdx(new).GetAtomicNum() <= 1:
            raise ValueError("Assembly port must lie on a heavy atom")
        counts[new] += 1
        ports.append({
            "parent_atom_index": old,
            "product_atom_index": new,
            "removed_parent_atom_index": edge["removed_atom"],
            "parent_bond_index": edge["bond_index"],
            "bond_order": 1,
        })
    tagged = Chem.Mol(product)
    for atom in tagged.GetAtoms():
        atom.SetAtomMapNum(2 ** counts[atom.GetIdx()] if atom.GetAtomicNum() > 1 else 0)
    text = Chem.MolToSmiles(Chem.RemoveHs(tagged), canonical=True, isomericSmiles=True)
    heavy = sum(a.GetAtomicNum() > 1 for a in product.GetAtoms())
    n = len(ports)
    return {
        "module_key": hashlib.sha256((MODULE_POLICY + "\0" + text).encode()).hexdigest(),
        "module_policy": MODULE_POLICY,
        "attachment_smiles": text,
        "heavy_atoms": heavy,
        "size_class": "small" if heavy <= 2 else "larger",
        "role": ("isolated" if n == 0 else "terminal" if n == 1
                 else "connector" if n == 2 else "junction"),
        "port_count": n,
        "port_atom_count": len(counts),
        "ports": ports,
    }
