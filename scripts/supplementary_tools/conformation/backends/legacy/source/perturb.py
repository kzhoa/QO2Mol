def bond_angle_extend_conformer(molecule,ignore_alkane=True):
    from .zmatrix import ZMatrix
    from ....graph import Element
    import math,copy

    _alkane_types = ["h_1", "c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4"]

    def is_atoms_alkane(m, *ids):
        return all(m.Atoms[i].atom_type_name in _alkane_types for i in ids)

    def is_in_same_ring(m, *ids):
        for k, v in m.ring_dict.items():
            if all(i in v[:-1] for i in ids):
                return True

        return False

    def determine_delta_bond(m, bond):
        if is_in_same_ring(m, bond.a1, bond.a2):
            return 0.05
        elif bond.type == "1":
            return 0.08
        elif bond.type == "3":
            return 0.05
        else:
            return 0.06

    def validate_pair_distance(mole, pos_array):
        for pair in getattr(mole, "Pair14", []) + getattr(mole, "Pair1n", []):
            elem1, elem2 = mole.Atoms[pair.a1].elem, mole.Atoms[pair.a2].elem
            vdw_diameter = Element.get(elem1).vdw_radius + Element.get(elem2).vdw_radius
            threshold = vdw_diameter / 2 ** (1 / 6)
            delta = pos_array[pair.a2] - pos_array[pair.a1]
            if delta.dot(delta) < threshold**2 and not is_in_same_ring(mole, pair.a1, pair.a2):
                return False

        return True
            

    zmat = ZMatrix()
    zmat.from_molecule(molecule)

    id_bonds = [
        i
        for i, b in enumerate(zmat.bonds)
        if b is not None and not (ignore_alkane and is_atoms_alkane(molecule, b.a1, b.a2))
    ]
    id_angles = [
        i
        for i, a in enumerate(zmat.angles)
        if a is not None
        and (a.type21 == "1" or a.type23 == "1")
        and not is_in_same_ring(molecule, a.a1, a.a3)
        and not (ignore_alkane and is_atoms_alkane(molecule, a.a1, a.a2, a.a3))
    ]

    _n_total = 0
    zmatrixes = {}

    half_len = math.ceil(len(id_bonds) / 2)
    for i, j in zip(id_bonds[:half_len], (id_bonds + [None])[half_len : 2 * half_len]):
        _n_total += 1
        zmatb = copy.deepcopy(zmat)
        zmatb.bonds[i].value += determine_delta_bond(molecule, zmatb.bonds[i])
        name = "b-%i" % (i + 1)
        if j is not None:
            zmatb.bonds[j].value -= determine_delta_bond(molecule, zmatb.bonds[j])
            name += "-%i" % (j + 1)
        if validate_pair_distance(molecule, zmatb.to_coordinates()):
            zmatrixes[name] = zmatb

    half_len = math.ceil(len(id_angles) / 2)
    for i, j in zip(id_angles[:half_len], (id_angles + [None])[half_len : 2 * half_len]):
        _n_total += 1
        zmata = copy.deepcopy(zmat)
        angle_i = zmata.angles[i]
        angle_i.value += 10 / 180 * math.pi
        angle_i.value_a += 10
        name = "a-%i" % (i + 1)
        if j is not None:
            angle_j = zmata.angles[j]
            angle_j.value -= 10 / 180 * math.pi
            angle_j.value_a -= 10
            name += "-%i" % (j + 1)
        if validate_pair_distance(molecule, zmata.to_coordinates()):
            zmatrixes[name] = zmata

    return zmatrixes

def ignore_alkane_torsion(molecule):
    scan_term = []
    for torsion in molecule.torsions:
        ats = [molecule.Atoms[an].atom_type_name for an in torsion]
        if not set(ats).issubset({"c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4", "h_1"}):
            scan_term.append(torsion)
    molecule.torsions = scan_term
