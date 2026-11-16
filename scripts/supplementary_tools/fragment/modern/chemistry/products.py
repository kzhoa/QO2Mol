"""Construct saturated boundaries in RDKit memory, without MOL text or files."""

from collections import Counter
from rdkit import Chem


def normalized_members(mol, components):
    raw = list(components)
    if not raw or any(
        type(i) is not int or not 0 <= i < mol.GetNumAtoms() for i in raw
    ):
        raise ValueError("Invalid fragment atom membership")
    members = set(raw)
    # A ring descriptor can include an exocyclic atom but omit its attached H.
    # Retain those original hydrogens (including isotope labels) before capping.
    members.update(
        n.GetIdx()
        for i in raw
        for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() == 1
    )
    return members


def connected_members(mol, members):
    remaining = set(members)
    parts = []
    while remaining:
        group, todo = set(), [min(remaining)]
        while todo:
            atom = todo.pop()
            if atom in group:
                continue
            group.add(atom)
            todo.extend(
                n.GetIdx()
                for n in mol.GetAtomWithIdx(atom).GetNeighbors()
                if n.GetIdx() in remaining and n.GetIdx() not in group
            )
        remaining -= group
        parts.append(sorted(group))
    return parts


def cap(mol, members, *, special_bond=()):
    members = set(members)
    if len(connected_members(mol, members)) != 1:
        raise ValueError("A fragment product must be connected")
    working = Chem.Mol(mol)
    Chem.Kekulize(working, clearAromaticFlags=True)
    deficits = Counter()
    boundary = []
    for bond in working.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a in members) != (b in members):
            order = int(bond.GetBondTypeAsDouble())
            if order not in (1, 2, 3):
                raise ValueError("Unsupported boundary bond order")
            kept, removed = (a, b) if a in members else (b, a)
            deficits[kept] += order
            boundary.append(
                dict(kept_atom=kept, removed_atom=removed, bond_order=order)
            )
    # The SSSR kernel identifies ring-internal bonds whose order increases when
    # another fused ring is removed. Account for that valence before adding H.
    special = set(special_bond) & members
    edits = []
    for bond in working.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in special and b in special:
            old = int(bond.GetBondTypeAsDouble())
            if old not in (1, 2) or deficits[a] < 1 or deficits[b] < 1:
                raise ValueError(
                    "SSSR bond-order adjustment exceeds available boundary valence"
                )
            bond.SetBondType(Chem.BondType.DOUBLE if old == 1 else Chem.BondType.TRIPLE)
            deficits[a] -= 1
            deficits[b] -= 1
            edits.append(dict(atoms=[a, b], old_order=old, new_order=old + 1))
    affected = {e["kept_atom"] for e in boundary}
    affected.update(i for e in edits for i in e["atoms"])
    cleared = []
    for atom in working.GetAtoms():
        atom.SetNoImplicit(True)
        if (
            atom.GetIdx() in affected
            and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        ):
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            cleared.append(atom.GetIdx())
    cleared_bonds = []
    for bond in working.GetBonds():
        if any(i not in members for i in bond.GetStereoAtoms()) or any(
            i in affected for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        ):
            if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                cleared_bonds.append(bond.GetIdx())
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            bond.SetBondDir(Chem.BondDir.NONE)
        if bond.GetBeginAtomIdx() in affected or bond.GetEndAtomIdx() in affected:
            bond.SetBondDir(Chem.BondDir.NONE)
    rw = Chem.RWMol(working)
    for old in reversed(range(mol.GetNumAtoms())):
        if old not in members:
            rw.RemoveAtom(old)
    ordered = sorted(members)
    mapping = {old: new for new, old in enumerate(ordered)}
    product_to_parent = list(ordered)
    caps = []
    for old in ordered:
        for _ in range(deficits[old]):
            hydrogen = Chem.Atom(1)
            hydrogen.SetNoImplicit(True)
            new = rw.AddAtom(hydrogen)
            rw.AddBond(mapping[old], new, Chem.BondType.SINGLE)
            product_to_parent.append(None)
            caps.append(dict(parent_atom_index=old, product_atom_index=new))
    product = rw.GetMol()
    product.RemoveAllConformers()
    Chem.SanitizeMol(product)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    # A successful SMILES export alone is insufficient: require atomwise
    # charge, isotope, radical and valence preservation across the boundary.
    for old, new in mapping.items():
        source, target = mol.GetAtomWithIdx(old), product.GetAtomWithIdx(new)
        before = (
            source.GetFormalCharge(),
            source.GetIsotope(),
            source.GetNumRadicalElectrons(),
            source.GetTotalValence(),
        )
        after = (
            target.GetFormalCharge(),
            target.GetIsotope(),
            target.GetNumRadicalElectrons(),
            target.GetTotalValence(),
        )
        if before != after:
            raise ValueError(
                f"Capping changed charge, isotope, radicals or valence at parent atom {old}: {before} -> {after}"
            )
    return product, dict(
        parent_atom_indices=ordered,
        product_to_parent=product_to_parent,
        boundary_bonds=boundary,
        hydrogen_caps=caps,
        caps=[[i, deficits[i], "H"] for i in ordered if deficits[i]],
        internal_bond_edits=edits,
        cleared_cut_stereocenters=cleared,
        cleared_double_bonds=cleared_bonds,
        coordinates_generated=False,
        normalization="rdkit_memory_boundary_v1",
        validation="atomwise_valence_charge_isotope_radical",
    )
