import itertools
import logging
from copy import deepcopy
import numpy as np
from ..q4_geometry import ME, MP
logger = logging.getLogger(__name__)

class ConformType:
    LOCAL_MINIMUM = "local minimum"
    STRETCH = "stretch"
    SCAN = "scan"
    CONSTRAINED_LOCAL_MINIMUM = "constrained local minimum"
    BARRIER = "barrier"
    ON_SLOPE = "on slope"
    OTHER = "other"
    SINGLE_POINT = "single point"
    OPTIMIZING = "optimizing"
    TORSION_SCAN_TYPES = [SCAN, CONSTRAINED_LOCAL_MINIMUM, BARRIER, ON_SLOPE, OTHER]

def local_minimum_pes(scan_curve):
    """
    找到一个分子local_minimum构象
    输入：
        scan_curve: List[Atom.No] 查找local_minimum的scan term
    输出：
        rlm_dict: Dict， 每个flexible torsion所包含的constrained local minimum值
    """
    def _base_molecule(scan_curve):

        rlm_dict = {}
        for m_name, terms in scan_curve.items():
            rlm_dict[m_name] = {}
            for term, molecules in terms.items():
                tmp = [molecule for molecule in molecules 
                       if molecule.conform_type in 
                       [ConformType.LOCAL_MINIMUM, ConformType.CONSTRAINED_LOCAL_MINIMUM]]

                tmp = sorted(tmp, key=lambda m:m.energy)
                rlm_dict[m_name][term] = [molecule.constrain[0].fix_value for molecule in tmp]
        return rlm_dict

    def _base_dict(scan_curve):
        rlm_dict = {}
        for m_name, terms in scan_curve.items():
            rlm_dict[m_name] = {}
            for term, data_dict_arr in terms.items():
                tmp = [data_dict for data_dict in data_dict_arr 
                       if data_dict["conform_type"] in 
                       [ConformType.LOCAL_MINIMUM, ConformType.CONSTRAINED_LOCAL_MINIMUM]]

                tmp = sorted(tmp, key=lambda m:m["energy"])
                rlm_dict[m_name][term] = [data_dict["constrain"][0][-1] for data_dict in tmp]
        return rlm_dict
    
    for m_name,vv in scan_curve.items():
        for torsion,vvv in vv.items():
            if isinstance(vvv[0],dict):
                return _base_dict(scan_curve)
            else:
                return _base_molecule(scan_curve)

def create_lm_by_combine_scan_rlm(molecule, rlm_dict_mol, n=64, create_constrain=False,idx=None):
    """
    输入：
        molecule: Molecule
        rlm_dict_mol: 从local_minimum_pes得到
        n: int, 生成local minimum结构不超过n
        create_constratin: 旋转的二面角是否变成constrain项
    输出：
        molecules: List[Molecule], 具有新坐标的molecule对象复制
    """
    molecule.update_topol_value()
    molecules, value_arr, rotation_atoms = [], [], []
    for aa, bb in rlm_dict_mol.items():
        value_arr.append(bb)
        atoms = [int(ss) for ss in aa.split("-")]
        rotation_atoms.append(atoms)
    for rr in itertools.islice(itertools.product(*value_arr), n):
        molecules.append(deepcopy(molecule))
        constrain_arr = []
        for i in range(len(rr)):
            if create_constrain:
                constrain_arr.append(
                    [
                        rotation_atoms[i][0],
                        rotation_atoms[i][1],
                        rotation_atoms[i][2],
                        rotation_atoms[i][3],
                        rr[i],
                    ]
                )

            ME._structure_change(molecules[-1],rotation_atoms[i],rr[i])
            
        if create_constrain:
            molecules[-1].create_constrain(constrain_arr)
    if idx is not None:
        return molecules, idx
    else:
        return molecules

def remove_similar_conformer(molecules,target_molecule=None):
    if target_molecule is not None:
        target_inertia = MP._inertia_calculate(target_molecule, ignore_hydrogen=True)[0]
    list_inertia = [MP._inertia_calculate(m, ignore_hydrogen=True)[0] for m in molecules]
    indx_keep = []

    for idx, inertia in enumerate(list_inertia):
        if target_molecule is not None:
            if np.allclose(inertia, target_inertia, rtol=0, atol=5.0):
                continue
        if any(np.allclose(inertia, list_inertia[ref], rtol=0, atol=5.0) for ref in indx_keep):
            continue
        indx_keep.append(idx)

    n_duplicate = len(molecules) - len(indx_keep)
    if n_duplicate:
        logger.debug(f"{n_duplicate}/{len(molecules)} duplicated conformations removed")

    return [molecules[ii] for ii in indx_keep]
