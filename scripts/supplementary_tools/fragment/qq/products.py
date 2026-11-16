"""Recompute real parent boundaries and replace removed single bonds with H."""

from rdkit import Chem


def cap(mol, members, *, allowed_ring_bonds=()):
    members = set(members)
    if not members or not members <= set(range(mol.GetNumAtoms())):
        raise ValueError("Invalid candidate atom membership")
    if any(type(i) is not int or not 0 <= i < mol.GetNumBonds() for i in allowed_ring_bonds):
        raise ValueError("Invalid authorized ring bond index")
    allowed_ring_bonds = set(allowed_ring_bonds)
    boundary = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (i in members) != (j in members):
            if (
                bond.GetBondType() != Chem.BondType.SINGLE
                or (bond.IsInRing() and bond.GetIdx() not in allowed_ring_bonds)
                or bond.GetIsAromatic()
            ):
                raise ValueError(
                    "Unsupported boundary: only non-ring, nonaromatic single bonds may be capped"
                )
            boundary.append(
                {
                    "bond_index": bond.GetIdx(),
                    "kept_atom": i if i in members else j,
                    "removed_atom": j if i in members else i,
                }
            )
    working = Chem.Mol(mol)
    cleared = []
    for edge in boundary:
        atom = working.GetAtomWithIdx(edge["kept_atom"])
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            cleared.append(atom.GetIdx())
    cleared_double_bonds = []
    for bond in working.GetBonds():
        if any(i not in members for i in bond.GetStereoAtoms()):
            if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                cleared_double_bonds.append(bond.GetIdx())
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            bond.SetBondDir(Chem.BondDir.NONE)
    if boundary:
        working = Chem.FragmentOnBonds(
            working, [b["bond_index"] for b in boundary], addDummies=True
        )
        for atom in working.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(1)
                atom.SetIsotope(0)
                atom.SetFormalCharge(0)
                atom.SetNumRadicalElectrons(0)
                atom.SetNoImplicit(True)
    mappings = []
    parts = Chem.GetMolFrags(
        working, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mappings
    )
    matches = [
        (part, mapping)
        for part, mapping in zip(parts, mappings)
        if {i for i in mapping if i < mol.GetNumAtoms()} == members
    ]
    if len(matches) != 1:
        raise ValueError("Candidate must be a single connected parent subgraph")
    product, mapping = matches[0]
    product.RemoveAllConformers()
    Chem.SanitizeMol(product)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    inverse = {old: new for new, old in enumerate(mapping)}
    caps = []
    for new, old in enumerate(mapping):
        if old >= mol.GetNumAtoms():
            kept = working.GetAtomWithIdx(old).GetNeighbors()[0].GetIdx()
            caps.append({"product_atom_index": new, "parent_atom_index": kept})
    if len(caps) != len(boundary):
        raise ValueError("Boundary/cap count mismatch")
    # Retain the parent's charge and radical state on every original atom.
    for old in members:
        a, b = mol.GetAtomWithIdx(old), product.GetAtomWithIdx(inverse[old])
        if (a.GetFormalCharge(), a.GetNumRadicalElectrons(), a.GetIsotope()) != (
            b.GetFormalCharge(),
            b.GetNumRadicalElectrons(),
            b.GetIsotope(),
        ):
            raise ValueError(
                "Capping changed input charge, radicals or isotope on a retained atom"
            )
    construction = {
        "parent_atom_indices": sorted(members),
        "boundary_bonds": boundary,
        "hydrogen_caps": caps,
        "product_to_parent": [i if i < mol.GetNumAtoms() else None for i in mapping],
        "cleared_cut_stereocenters": sorted(set(cleared)),
        "cleared_double_bonds": cleared_double_bonds,
        "coordinates_generated": False,
    }
    opened = [b['bond_index'] for b in boundary if mol.GetBondWithIdx(b['bond_index']).IsInRing()]
    if opened:
        construction['opened_ring_bonds'] = opened
        construction['ring_membership_changed_parent_atoms'] = sorted(
            old for old in members if mol.GetAtomWithIdx(old).IsInRing()
            != product.GetAtomWithIdx(inverse[old]).IsInRing())
    return product, construction
