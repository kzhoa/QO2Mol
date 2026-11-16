def flexible_torsion(molecule,idx=None):
    ignore_atom_style = ["dummy"]
    if len(molecule.Atoms) < 2:
        return None
    for bb in molecule.Bonds:
        a1 = bb.a1
        a2 = bb.a2
        if (
            molecule.Atoms[a1].bond_type[molecule.Atoms[a1].connect.index(a2)] == "1"
            and molecule.Atoms[a1].s not in ignore_atom_style
            and molecule.Atoms[a2].s not in ignore_atom_style
        ):
            if molecule.Atoms[a1].connectivity_type[molecule.Atoms[a1].connect.index(a2)] in [
                "S",
                "eS",
                "bS",
                "br",
                "J",
                "eJ",
                "bJ",
            ]:
                if molecule.Atoms[a1].local != "EN" and molecule.Atoms[a2].local != "EN":
                    if molecule.Atoms[a1].local not in ["LT", "LTC"] and molecule.Atoms[a2].local not in ["LT", "LTC"]:
                        bb.flexible = "yes"
                    else:
                        if len(molecule.Atoms[a1].connect) > 1 and len(molecule.Atoms[a2].connect) > 1:
                            cea = None
                            if molecule.Atoms[a1].local in ["LT", "LTC"]:
                                cea = a1
                            else:
                                if molecule.Atoms[a2].local in ["LT", "LTC"]:
                                    cea = a2
                            if cea is not None:
                                if molecule.Atoms[cea].elem != "C":
                                    bb.flexible = "yes"
                                else:
                                    elems = []
                                    for ca in molecule.Atoms[cea].connect:
                                        if ca != a1 and ca != a2:
                                            elems.append(molecule.Atoms[ca].elem)
                                    if len(set(elems)) > 1:
                                        bb.flexible = "yes"   
    if idx is not None:     
        return molecule,idx

def scan_torsion(molecule, ignore_alkane=False,idx=None):
    """
    确定scan_torsion。
    输入：
        m: Molecule
        ignore_alkane: True or False。是否忽略烷基链部分的二面角
    输出：
        scan_term: List[[int]]可以被旋转的二面角
    """
    ignore_atom_style = ["dummy"]
    scan_term = []
    if len(molecule.Atoms) < 4:
        molecule.torsions = []
    for bb in molecule.Bonds:
        if not hasattr(bb, "flexible") or bb.flexible != "yes":
            continue
        a1 = bb.a1
        a2 = bb.a2
        if "3" in molecule.Atoms[a1].bond_type or "3" in molecule.Atoms[a2].bond_type:
            continue
        #ats = [molecule.Atoms[a1].atom_type_name, molecule.Atoms[a2].atom_type_name]
        neighbors1 = [molecule.Atoms[i] for i in molecule.Atoms[a1].connect if i != a2 and molecule.Atoms[i].s not in ignore_atom_style]
        atom0 = neighbors1[0]
        for a in neighbors1[1:]:
            if a.mass > atom0.mass:
                atom0 = a
        #ats.append(atom0.atom_type_name)
        neighbors2 = [molecule.Atoms[i] for i in molecule.Atoms[a2].connect if i != a1 and molecule.Atoms[i].s not in ignore_atom_style]
        atom3 = neighbors2[0]
        for a in neighbors2[1:]:
            if a.mass > atom3.mass:
                atom3 = a
        #ats.append(atom3.atom_type_name)
        #if ignore_alkane and set(ats).issubset({"c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4", "h_1"}):
        #    continue
        scan_term.append([atom0.No, a1, a2, atom3.No])
    molecule.torsions = scan_term
    if idx is not None:
        return molecule,idx
