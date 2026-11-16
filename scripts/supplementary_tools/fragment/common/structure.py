# Shared structural helper extracted from the original implementation.
from copy import deepcopy
from .constants import *
SCAOG = SINGLE_CONNECT_ATOM_OR_GROUP

def get_2d_connectivity(molecule):
    element_number = {ii:elem for ii,elem in enumerate(molecule.elements)}
    all_connect = {ii:connect for ii,connect in enumerate(molecule.connectivity)}
    return element_number, all_connect

def remove_hydrogen_atoms(element_number,all_connect):
    reduce_connect = deepcopy(all_connect)
    H_atom_arr = []
    for item in all_connect.keys():
        if element_number[item] in SCAOG:
            del reduce_connect[item]
            H_atom_arr.append(item)
    for item in reduce_connect.keys():
        tmp = reduce_connect[item][:]
        n = len(reduce_connect[item])
        for i in range(0, n):
            if reduce_connect[item][i] in H_atom_arr:
                tmp.remove(reduce_connect[item][i])
        reduce_connect[item] = tmp[:]
    terminal_atom_arr = []
    for item in reduce_connect.keys():
        if len(reduce_connect[item]) == 1:
            terminal_atom_arr.append(int(item))
    return reduce_connect, terminal_atom_arr, H_atom_arr

def assign_structure_info(molecule, rings,ring_blocks,ring_block_components):
    #rings = {}
    #ring_stru = []
    #for arr in cyclos:
    #    n = len(arr) - 1
    #    nn = int(n / 2)
    #    ss = f"R{n}-{arr[0]}_{arr[nn]}_{arr[-1]}"
    #    if ss in rings.keys():
    #        ss = f"R{n}-{arr[0]}_{arr[nn]}--{randint(0, 100)}_{arr[-1]}"
    #    rings[ss] = arr
#
    #tmp = []
    #for aa, bb in rings.items():
    #    tmp.append(bb[:-1])
    #    ring_stru = combine_arr(tmp)

    for i in range(0, len(molecule.Atoms)):
        inring = []
        ring_size = []
        ring_prop = []
        for aa, bb in rings.items():
            if i in bb:
                inring.append(aa)
                ring_size.append(len(bb) - 1)
                ring_prop.append(bb[-1])
        molecule.Atoms[i].has_ring = inring
        molecule.Atoms[i].has_ring_size = ring_size
        molecule.Atoms[i].has_ring_property = ring_prop

    molecule.rings = rings
    molecule.ring_blocks = ring_blocks
    molecule.ring_block_components = ring_block_components

def assign_local_info(molecule,reduce_connect):
    cycloset = []
    for __, cyclo in molecule.ring_dict.items():
        cycloset.extend(cyclo[:-1])

    cycloset = set(cycloset)

    for atom_n,atom in enumerate(molecule.Atoms):
        if atom.elem in SCAOG:
            atom.local = "EN"
            continue

        connect_atom_number_in_ring = len(set(reduce_connect[atom_n]) & cycloset)
        connect_atom_number_out_ring = len(reduce_connect[atom_n]) - connect_atom_number_in_ring
        
        if atom_n not in cycloset:
            if len(reduce_connect[atom_n]) == 1:
                if connect_atom_number_in_ring == 0:
                    atom.local = "LT"
                else:
                    atom.local = "LTC"
            else:
                if connect_atom_number_in_ring == 0:
                    atom.local = "LM"
                else:
                    atom.local = "LBC" if connect_atom_number_in_ring > 1 else "LMC"
            continue

        if len(atom.inring) == 1:
            if connect_atom_number_in_ring == 2:
                atom.local = "CL" if connect_atom_number_out_ring != 0 else "C"
            else:
                atom.local = "CBL" if connect_atom_number_out_ring != 0 else "CB"
            continue

        bridge = 0
        for connect_n in atom.connect:
            if connect_n in cycloset and len(set(molecule.Atoms[connect_n].inring) & set(atom.inring)) == 0:
                bridge += 1
        if connect_atom_number_in_ring - bridge == 2:
            atom.local = "CFM" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CFMB" if bridge != 0 else "CFML")
        elif connect_atom_number_in_ring - bridge == 3:
            atom.local = "CF" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CFB" if bridge != 0 else "CFL")
        else:
            atom.local = "CS" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CSB" if bridge != 0 else "CSL")

def get_atom_flag(x):
    if x in LOCAL_LINK:
        return "c"
    else:
        if x in LOCAL_CYCLE_F:
            return "f"
        else:
            if x in LOCAL_CYCLE_B:
                return "b"
            else:
                return "x"

def get_connect(bn,str,c1,c2):
    if str == "bb":
        str = "nr" if len(set(c1).intersection(c2)) == 0 else "rr"
    index_label = [["cc"],["cb","cf","cx","bc","fc","xc"],["ff"],["nr"],
             ["xx","bf","bx","fb","fx","xf","xb","rr"]]
    for i,arr in enumerate(index_label):
        if str in arr:
            index = i + 1
            break
    types = {
            1: ["S", "D", "T", "J", "M"], #1
            2: ["eS", "eD", "eM", "es", "eM"],#2
            3: ["fs", "fd", "ft", "fr", "fm"],#3
            4: ["bS", "bD", "bm", "br", "bm"],#4
            5: ["s", "d", "t", "r", "m"],#5
            }
    btype = ["1", "2", "3", "ar"]
    if bn not in btype:
        return types[index][-1]
    else:
        return types[index][btype.index(bn)]

def get_aromatic_connect(molecule):
    arom_atom_dict = {}
    for __, ring in molecule.rings.items():
        if ring[-1] in ["ar1", "ar2"]:
            for ai in ring[:-1]:
                if ai not in arom_atom_dict:
                    arom_atom_dict[ai] = []
                arom_atom_dict[ai] = list(set(arom_atom_dict[ai]).union(set(ring[:-1])))

    for atom in molecule.Atoms:
        #if not hasattr(atom, "bond_type_old"):
        #    atom.bond_type_old = deepcopy(atom.bond_type)
        if atom.ID in arom_atom_dict:
            atom.bond_type_aromatic = []
            for i in range(len(atom.connectivity)):
                if atom.connectivity[i] in arom_atom_dict[atom.No]:
                    atom.bond_type_aromatic.append("ar")
                else:
                    atom.bond_type_aromatic.append(atom.bond_type[i])
        else:
            atom.bond_type_aromatic = deepcopy(atom.bond_type)

def assign_connect_type(molecule):
    get_aromatic_connect(molecule)
    for atom_n, atom in enumerate(molecule.Atoms):
        atom_cyclo = atom.inring
        atom_str = get_atom_flag(atom.local)
        atom.connectivity_type = []
        for atom_m, conn_an in enumerate(atom.connect):
            atom_conn = molecule.Atoms[conn_an]
            atom_conn_cyclo = atom_conn.inring
            atom_conn_str = get_atom_flag(atom_conn.local)
            bn = atom.bond_type_aromatic[atom_m]
            atom.connectivity_type.append(get_connect(bn,f"{atom_str}{atom_conn_str}",atom_cyclo,atom_conn_cyclo))

def assign_conjugate_info(molecule):
    for atom in molecule.Atoms:
        atom.bond_type_conjugate = ["" for _ in atom.connect]
    conju_atoms = []
    for ii,atom in enumerate(molecule.Atoms):
        conju_flag = False
        if set(atom.connectivity_type) & set(MULTI_BOND) == set():
            continue
        for jj,cbn in enumerate(atom.connectivity_type):
            atom1 = molecule.Atoms[atom.connectivity[jj]]
            if cbn in SINGLE_BOND and len(set(atom1.connectivity_type) & set(MULTI_BOND)) > 0:
                conju_flag = True
            if conju_flag:
                break
        
        if conju_flag:
            conju_atoms.append(ii)

    for ii,atom in enumerate(molecule.Atoms):
        for jj,an in enumerate(atom.connectivity):
            if an in conju_atoms:
                if atom.connectivity_type[jj] in CHAIN_BOND:
                    atom.bond_type_conjugate[jj] = "J"
                elif atom.connectivity_type[jj] in RING_BOND:
                    atom.bond_type_conjugate[jj] == "j"

def divide_charge(molecule):
    """
    Sum out formal charges between bonded opposite charges
    Average formal charges between several uni-coordinated O/S on C/S/P
    Average formal charges between several tri-coordinated N on C
    Spread formal charges to all heteroatoms in a conjugated system (TODO Only N-containing 5-member aromatic rings are considered)
    """
    charge_group = []
    for atom in molecule.Atoms:
        atom.primitive_formal_charge = atom.formal_charge
    for bond in molecule.Bonds:
        atom1, atom2 = molecule.Atoms[bond.a1], molecule.Atoms[bond.a2]
        if {atom1.formal_charge, atom2.formal_charge} == {1, -1}:
            atom1.primitive_formal_charge = 0
            atom2.primitive_formal_charge = 0
    _processed = []
    for atom in molecule.Atoms:
        if atom in _processed:
            continue
        if atom.elem in ["O", "S"] and len(atom.connect) == 1 and atom.primitive_formal_charge == -1:
            root = molecule.Atoms[atom.connect[0]]
            if root.elem not in ["C", "S", "P"]:
                continue
            neighbors = [
                molecule.Atoms[i]
                for i in root.connect
                if molecule.Atoms[i].elem in ["O", "S"] and len(molecule.connectivity[i]) == 1
            ]
            if len(neighbors) > 1:
                charge = sum(neigh.primitive_formal_charge for neigh in neighbors)
                for neigh in neighbors:
                    neigh.primitive_formal_charge = charge / len(neighbors)
                _processed.extend(neighbors)
                if charge > 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID for atom in neighbors]+[root.ID,label])

    #####改该函数必测试：[O-]CC(N1N=N[NH+]=C(N)1)C([O-])=O,result:[10, 11, 9], [5, 4, 3, 2, 0, 1], [8]
    rn_used = []
    _processed = []
    for atom in molecule.Atoms:
        if atom in _processed:
            continue
        if atom.elem in ["N"] and len(atom.connect) == 3 and atom.primitive_formal_charge == 1:
            for idx_neigh in atom.connect:
                root = molecule.Atoms[idx_neigh]
                if root.elem not in ["C"]:
                    if atom.formal_charge> 0:
                        label = "positive"
                    else:
                        label = "negative"
                    charge_group.append([atom.ID,label])
                    continue
                neighbors = [
                    molecule.Atoms[i]
                    for i in root.connect
                    if molecule.Atoms[i].elem in ["N"] and len(molecule.connectivity[i]) == 3
                ]
                if len(neighbors) > 1:
                    for ring in atom.has_ring:
                        ss = ring.split("_")
                        if ring[-3:] in ["ar1","ar2"] and ring[1] in ["5","6"]:
                            members = [molecule.Atoms[an] for an in molecule.rings[ring][:-1]]
                            neighbors += [atom for atom in members if atom.elem == "N"]
                            rn_used.append(ring)
                    neighbors = list(set(neighbors))
                    charge = sum(neigh.primitive_formal_charge for neigh in neighbors)
                    for neigh in neighbors:
                        neigh.primitive_formal_charge = charge / len(neighbors)
                    _processed.extend(neighbors)
                    if charge > 0:
                        label = "positive"
                    else:
                        label = "negative"
                    charge_group.append([atom.ID for atom in neighbors]+[root.ID,label])
                    break
    
    for rn,ring in molecule.ring_dict.items():
        if len(ring) != 6 or ring[-1] not in ["ar1", "ar2"]:
            continue
        if rn not in rn_used:
            members = [molecule.Atoms[i] for i in ring[:-1]]
            nitrogens = [atom for atom in members if atom.elem == "N"]
            charge = sum(atom.primitive_formal_charge for atom in nitrogens)
            if charge != 0:
                if charge > 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID for atom in nitrogens])
            for atom in nitrogens:
                atom.primitive_formal_charge = charge / len(nitrogens)

    _tmp = []
    for rr in charge_group:
        _tmp.extend(rr)
    for atom in molecule.Atoms:
        if atom.formal_charge != 0:
            if atom.ID not in _tmp:
                if  atom.formal_charge> 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID,label])

    molecule.charge_group = charge_group

def assign_hybrid(molecule,idx=None):
    _label = {
            "C":{0:"s",1:"sp",2:"sp",3:"sp2",4:"sp3"},
            "H":{0:"s",1:"s",},
            "O":{0:"s",1:"sp2",2:"sp3"},
            "N":{0:"s",1:"sp",2:"sp2",3:"sp3",4:"sp3+",5:"sp3d2"},
            "S":{0:"s",1:"sp2",2:"sp3",3:"sp2d",4:"sp3d2"},
            "P":{0:"s",1:"sp",2:"sp2",3:"sp3",4:"sp3+",5:"sp3d2"},
            "F":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Cl":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q"},
            "Br":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q"},
            "Si":{0:"s",1:"s",2:"sp",3:"sp2",4:"sp3"},
            "B":{0:"s",1:"s",2:"sp2",3:"sp3",4:"sp3-"},
            "I":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q",5:"sp3f"},
            "Na":{-1:"s",0:"s",1:"s",2:"s",3:"s",4:"s"},
            "K":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Ca":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Mg":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            }
    for atom in molecule.Atoms:
        bond_number = len(atom.connectivity)
        bond_number -= atom.formal_charge

        if atom.element in ["N","P"]:
            if len(atom.connectivity) == 4:
                bond_number = 4
            else:
                o_num = sum([int(bn) - molecule.Atoms[atom.connectivity[jj]].formal_charge 
                             for jj,bn in enumerate(atom.bond_type) 
                             if molecule.Atoms[atom.connectivity[jj]].element == "O"])
                if o_num >= 2:
                    bond_number = 5
                    #bond_number += o_num 
                    #if bond_number >= 4:
                    #    bond_number = 5
        if atom.element == "B" and len(atom.connectivity) == 4 and atom.formal_charge == -1:
            # Four-coordinate borate is explicitly covered by the typing
            # resource; subtracting its charge must not create a fifth ligand.
            atom.hybrid = "sp3-"
        else:
            try:
                atom.hybrid = _label[atom.element][bond_number]
            except KeyError as exc:
                raise ValueError(
                    f"Legacy hybridization table does not support atom {atom.ID} "
                    f"({atom.element}, coordination={len(atom.connectivity)}, charge={atom.formal_charge})"
                ) from exc
    if idx is not None:
        return molecule,idx
    else:
        return molecule
