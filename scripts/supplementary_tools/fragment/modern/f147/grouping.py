"""147 grouping rules specialized for fragmentation."""

from copy import deepcopy
from pathlib import Path
from .typing import AtomType
from .connectivity import find_mole
from ..chemistry.chain import search_chain
from ..chemistry.elements import Element
from ..chemistry.constants import (
    SINGLE_CONNECT_ATOM_OR_GROUP as SCAOG,
    BREAK_BOND,
    LOCAL_CYCLE,
)

LOCAL_EN = ["EN", "LT", "LTC"]
atf = str(Path(__file__).with_name("fg_def.json"))


class Grouping:

    def __init__(self):
        self._init_setting()

    def _init_setting(self):
        self._fgs = {
            "inorganic": [
                "HCl",
                "HF",
                "HBr",
                "HI",
                "Cl-",
                "F-",
                "Br-",
                "I-",
                "Na+",
                "K+",
                "Ca2+",
                "Mg2+",
            ],
            "rings": [],
            "bsi": ["b4neu", "b4XN", "b4FN", "SiX", "SiF", "BO", "B", "SiO", "Si"],
            "ion1": ["p4p", "p3p", "amo", "hamp", "amp", "imp"],
            "cum": [
                "numno",
                "numns",
                "numnn",
                "numo",
                "nums",
                "numn",
                "num",
                "cumno",
                "cumns",
                "cumnn",
                "cumo",
                "cums",
                "cumn",
                "cum",
            ],
            "tribond": ["da", "ic", "cy"],
            "az": ["az"],
            "hacid": [
                "hna",
                "hnf",
                "hni",
                "hsa",
                "hsai",
                "hsah",
                "hsf",
                "hsfi",
                "hsfh",
                "hso",
                "hsi",
                "hsif",
                "hsio",
                "hpa",
                "hpai",
                "hpah",
                "hpf",
                "hpfi",
                "hpfh",
                "hpo",
                "hpoi",
                "hpoh",
                "hpi",
                "hp3a",
                "hpif",
                "hpio",
            ],
            "acid": [
                "na",
                "nf",
                "ni",
                "nif",
                "n3a",
                "n3f",
                "sa",
                "sai",
                "sah",
                "sf",
                "sfi",
                "sfh",
                "so",
                "si",
                "sif",
                "sio",
                "pa",
                "pai",
                "pah",
                "pf",
                "pfi",
                "pfh",
                "po",
                "poi",
                "poh",
                "pi",
                "p3a",
                "pif",
                "pio",
            ],
            "ion2": ["ai", "ti"],
            "b0": [
                "_ox_",
                "_hz_",
                "_ke_",
                "_im_",
                "_tk_",
                "ox_",
                "hz_",
                "ke_",
                "im_",
                "tk_",
            ],
            "b": ["ox", "hz", "is", "ke", "im", "tk"],
            "a": ["al", "am", "th", "ph"],
            "a1": ["et", "te", "sm", "tm"],
            "halogen": ["f", "cl", "br", "i"],
            "yne": ["yne_conju", "yne"],
            "ene": ["ene_conju", "ene"],
            "ane": ["ane"],
            "res": [],
            "double": ["ene", "yne", "ene_conju", "yne_conju"],
            "ch_conju": ["yne_conju", "ene_conju"],
        }
        self._create_res_fgs()
        self._fgs["a_type"] = (
            self._fgs["ion2"] + self._fgs["a"] + self._fgs["a1"] + self._fgs["halogen"]
        )
        self._fgs["couple_type"] = []
        attrs = ["tribond", "b", "b0", "yne", "ene", "az"]
        for attr in attrs:
            self._fgs["couple_type"].extend(self._fgs[attr])
        self._fgs["double_type"] = self._fgs["ene"] + self._fgs["yne"] + self._fgs["az"]
        self._fgs["type_0"] = self._fgs["inorganic"] + self._fgs["bsi"]
        self._fgs["type_1"] = self._fgs["acid"] + self._fgs["hacid"]
        self._fgs["type_2"] = self._fgs["cum"]
        self._fgs["type_3"] = self._fgs["ion1"]
        self._fgs["type_4"] = self._fgs["b0"] + self._fgs["b"]
        self._fgs["type_5"] = (
            self._fgs["tribond"] + self._fgs["az"] + self._fgs["yne"] + self._fgs["ene"]
        )
        self._fgs["type_6"] = []
        self._fgs["type_7"] = self._fgs["a_type"]
        self._fgs["type_8"] = self._fgs["ane"]
        self._a_label = {
            "al": "o",
            "ai": "o",
            "et": "o",
            "th": "s",
            "te": "s",
            "ti": "s",
            "am": "n",
            "sm": "n",
            "tm": "n",
            "ph": "p",
            "f": "f",
            "cl": "cl",
            "br": "br",
            "i": "i",
            "o": "o",
            "s": "s",
            "n": "n",
            "p": "p",
        }
        self._a_symbol = {
            "ai": "a",
            "al": "b",
            "et": "c",
            "am": "d",
            "sm": "e",
            "tm": "f",
            "ti": "g",
            "th": "i",
            "te": "j",
            "ph": "k",
            "o": "a",
            "n": "b",
            "s": "c",
            "p": "d",
            "O": "a",
            "N": "b",
            "S": "c",
            "P": "d",
            "f": "l",
            "cl": "m",
            "br": "n",
            "i": "o",
            "F": "l",
            "Cl": "m",
            "Br": "n",
            "I": "o",
            "": "z",
        }
        self._tag_label = {
            1: "alpha",
            2: "beta",
            3: "gamma",
            4: "delta",
            5: "epsilon",
            6: "zeta",
            7: "eta",
            8: "theta",
            9: "lota",
            10: "kappa",
            11: "lambda",
            12: "mu",
            13: "nu",
            14: "xi",
            15: "omicron",
            16: "pi",
            17: "rho",
            18: "sigma",
            19: "tau",
            20: "upsilon",
            21: "phi",
            22: "chi",
            23: "psi",
            24: "omega",
        }
        self._fg_order = []
        for i in range(0, 9):
            self._fg_order.extend(self._fgs["type_%d" % i])
        self.at = AtomType(atf=atf, check_at_types=False)

    def _create_res_fgs(self):
        __elems = [
            "ONSP",
            "NSP",
            "OSP",
            "ONP",
            "ONS",
            "SP",
            "NP",
            "NS",
            "OP",
            "OS",
            "ON",
            "O",
            "N",
            "S",
            "P",
            "",
        ]
        __ars = ["ar1", "ar2", "ar3", "ar4", "ar5", ""]
        __rings = ["3R", "4R", "5R", "6R", "7R", "8R", "9R"]
        ring_fg = []
        for rr in __rings:
            for elems in __elems:
                for ar in __ars:
                    if ar != "":
                        ring_fg.append(f"{rr}{elems}_{ar}")
                    else:
                        ring_fg.append(f"{rr}{elems}")
        self._fgs["rings"] = ring_fg
        self._fgs["rings"].append("MR")
        self._fgs["rings"].append("LR")
        self._fgs["rings"].append("SR")
        self._fgs["rings"].append("R")

    def _clean(self, molecule):
        for atom in molecule.Atoms:
            for attr in [
                "function_group",
                "function_group_id",
                "function_group_tag",
                "_function_group_reverse",
                "_combine_a_fg_order",
                "_a_b_a_fg_order",
                "_function_group_direction",
                "_double_fg_order",
            ]:
                if hasattr(atom, attr):
                    delattr(atom, attr)
        for attr in ["_anhy_fgs", "_double_fgs", "_a_b_a_fgs", "_combine_a_fgs"]:
            if hasattr(molecule, attr):
                delattr(molecule, attr)

    def _assign_atom_type(self, molecule, atom_arr=None):
        self.at._atom_type(molecule, atoms_arr=atom_arr)

    def _assign_connect_fg(self, fg, fgn, bt, atom, post):
        _tag = {"1": "s", "2": "d", "3": "t"}
        atom.function_group_tag = f"{atom.elem}_{_tag[bt]}_{post}"
        atom.function_group = fg
        atom.function_group_id = fgn

    def _get_chain_center(self, molecule, res_atoms):
        fg_center = []
        not_a_center = []
        a_center_dict = {}
        couple_pair = {}
        _tmp = []
        not_center = []
        for atom in molecule.Atoms:
            if atom.ID in res_atoms and ":" in atom.atom_type_name:
                this_fg, ee = atom.atom_type_name.split(":")
                flag = True
                if this_fg in self._fgs["couple_type"]:
                    if atom.ID not in _tmp:
                        couple_pair[atom.ID] = [
                            bb
                            for jj, bb in enumerate(atom.connectivity)
                            if atom.bond_type[jj] != "1"
                        ][0]
                        _tmp.append(atom.ID)
                        _tmp.append(couple_pair[atom.ID])
                    else:
                        flag = False
                if flag:
                    this_index = self._fg_order.index(this_fg)
                    fg_center.append([atom.ID, this_index])
                    _arr_tmp = self._fgs["a_type"]
                    if this_fg not in _arr_tmp:
                        not_a_center.append(atom.ID)
                    else:
                        a_center_dict[atom.ID] = this_fg
            else:
                not_center.append(atom.ID)
        fg_center = sorted(fg_center, key=lambda x: x[1])
        return (fg_center, not_a_center, couple_pair, a_center_dict)

    def _assign_sing_chain_fg_acid_type(
        self, molecule, a, fgn, fg, tag, arr, not_a_center
    ):
        _tmp = []
        arr.append(a)
        _tmp.append(a)
        _tmp_atoms = []
        molecule.Atoms[a].function_group = fg
        molecule.Atoms[a].function_group_tag = f"{tag}_center"
        molecule.Atoms[a].function_group_id = fgn
        for ii, an in enumerate(molecule.Atoms[a].connectivity):
            if an not in arr:
                atom = molecule.Atoms[an]
                if molecule.Atoms[a].bond_type[ii] == "1":
                    if an not in not_a_center:
                        arr.append(an)
                        _tmp.append(an)
                        _tmp_atoms.append(atom)
                        atom.function_group_tag = f"{atom.element}_center"
                        self._assign_connect_fg(
                            fg, fgn, molecule.Atoms[a].bond_type[ii], atom, "alpha"
                        )
                elif atom.atom_type_name.split(":")[0] not in self._fgs["type_1"]:
                    arr.append(an)
                    _tmp.append(an)
                    self._assign_connect_fg(
                        fg, fgn, molecule.Atoms[a].bond_type[ii], atom, "alpha"
                    )
                    for jj, an1 in enumerate(atom.connectivity):
                        atom1 = molecule.Atoms[an1]
                        if an1 not in arr and ":" not in atom1.atom_type_name:
                            arr.append(an1)
                            _tmp.append(an1)
                            self._assign_connect_fg(
                                fg, fgn, atom.bond_type[jj], atom1, "beta"
                            )
        _o_type = [
            atom.atom_type_name.split(":")[0]
            for atom in _tmp_atoms
            if hasattr(atom, "atom_type_name")
        ]
        if fg in ["sa", "sf", "hsa", "hsf", "pa", "hpa", "pf", "hpf", "po", "hpo"]:
            if "ai" in _o_type or "ti" in _o_type:
                fg = f"{fg}i"
            elif "al" in _o_type:
                fg = f"{fg}h"
        for an in _tmp:
            molecule.Atoms[an].function_group = fg
        return arr

    def _assign_sing_chain_fg_a_type(self, molecule, a, fgn, fg, tag, arr):
        arr.append(a)
        molecule.Atoms[a].function_group = fg
        molecule.Atoms[a].function_group_tag = f"{tag}_center"
        molecule.Atoms[a].function_group_id = fgn
        return arr

    def _assign_sing_chain_fg_non_extend_type(
        self, molecule, a, fgn, fg, tag, arr, not_a_center
    ):
        arr.append(a)
        molecule.Atoms[a].function_group = fg
        molecule.Atoms[a].function_group_tag = f"{tag}_center"
        molecule.Atoms[a].function_group_id = fgn
        for ii, an in enumerate(molecule.Atoms[a].connectivity):
            if an not in arr:
                atom = molecule.Atoms[an]
                if an not in not_a_center:
                    arr.append(an)
                    self._assign_connect_fg(
                        fg, fgn, molecule.Atoms[a].bond_type[ii], atom, "alpha"
                    )
        return arr

    def _assign_single_chain_fg_b_extend_type(
        self, molecule, a, fgn, fg, tag, arr, couple_pair, not_a_center, a_center_dict
    ):

        def _get_fg_order(_extend, fg):
            _extend_tmp = [
                [[[], "", "z"], [[], "", "z"]],
                [[[], "", "z"], [[], "", "z"]],
                [[[], "", "z"], [[], "", "z"]],
                [[[], "", "z"], [[], "", "z"]],
                [[[], "", "z"], [[], "", "z"]],
            ]
            for ii, rr in enumerate(_extend):
                for jj, rrr in enumerate(rr):
                    _extend_tmp[ii][jj] = [[rrr[0]], rrr[1], self._a_symbol[rrr[1]]]
            _extend_ = []
            for ii, rr in enumerate(_extend_tmp):
                _extend_.append(
                    [
                        rr[0][0] + rr[1][0],
                        f"{rr[0][1]}{rr[1][1]}",
                        f"{rr[1][1]}{rr[0][1]}",
                    ]
                )
                if rr[1][2] != "z":
                    _extend_[-1].append(f"{rr[0][2]}{rr[1][2]}")
                else:
                    _extend_[-1].append(f"{rr[1][2]}{rr[0][2]}")
            _extend_ = sorted(_extend_, key=lambda x: x[-1])
            if len(_extend_[1][0]) == 0:
                if len(_extend_[0][0]) == 0:
                    _fg = fg
                    _fg_reverse = None
                    _fg_order_index = 5
                else:
                    _fg = f"{fg}_{_extend_[0][1]}"
                    _fg_reverse = f"{_extend_[0][2]}_{fg}"
                    if "z" in _extend_[0][3]:
                        _fg_order_index = 4
                    else:
                        _fg_order_index = 3
                    _extend_[0].append("right")
            else:
                _fg = f"{_extend_[0][2]}_{fg}_{_extend_[1][1]}"
                _fg_reverse = f"{_extend_[1][2]}_{fg}_{_extend_[0][1]}"
                if "z" in _extend_[0][3]:
                    _fg_order_index = 2
                elif "z" in _extend_[1][3]:
                    _fg_order_index = 1
                else:
                    _fg_order_index = 0
                _extend_[0].append("left")
                _extend_[1].append("right")
            _fg_order = f"{_fg_order_index}:{self._fg_order.index(fg)}:{_extend_[0][-1]}{_extend_[1][3]}"
            return ([_extend_[0], _extend_[1]], _fg, _fg_reverse, _fg_order)

        arr.append(a)
        molecule.Atoms[a].function_group_tag = f"{tag}_center"
        catom = molecule.Atoms[couple_pair[a]]
        _tmp = [a]
        _tmp.append(couple_pair[a])
        catom.function_group_tag = f"{catom.element}_fgcore"
        arr.append(couple_pair[a])
        if fg == "im":
            for ii, an in enumerate(catom.connectivity):
                if an not in arr and an not in not_a_center:
                    if molecule.Atoms[an].element == "O":
                        fg = "ox"
                    elif molecule.Atoms[an].element == "N":
                        fg = "hz"
                    elif molecule.Atoms[an].element == "S":
                        fg = "is"
                    _tmp.append(an)
                    arr.append(an)
                    molecule.Atoms[an].function_group_tag = (
                        f"{molecule.Atoms[an].element}_fgcore_extend"
                    )
                    break
        _extend = [[] for _ in range(len(molecule.Atoms[a].connectivity))]
        for ii, an in enumerate(molecule.Atoms[a].connectivity):
            if an not in arr and an not in not_a_center and (an in a_center_dict):
                atom = molecule.Atoms[an]
                _tmp.append(an)
                arr.append(an)
                atom.function_group_tag = f"{atom.element}_alpha"
                _extend[ii].append([an, a_center_dict[an]])
                for jj, an1 in enumerate(atom.connectivity):
                    if (
                        an1 not in arr
                        and an1 not in not_a_center
                        and (an1 in a_center_dict)
                    ):
                        atom1 = molecule.Atoms[an1]
                        _tmp.append(an1)
                        arr.append(an1)
                        atom1.function_group_tag = f"{atom1.element}_beta"
                        _extend[ii].append([an1, a_center_dict[an1]])
        _extend = [rr for rr in _extend if len(rr) > 0]
        _extend_, fg, fg_reverse, fg_order = _get_fg_order(_extend, fg)
        for an in _tmp:
            molecule.Atoms[an].function_group = fg
            molecule.Atoms[an].function_group_id = fgn
            if fg_reverse is not None:
                molecule.Atoms[an]._function_group_reverse = fg_reverse
                molecule.Atoms[an]._a_b_a_fg_order = fg_order
        for rr in _extend_:
            for an in rr[0]:
                molecule.Atoms[an]._function_group_direction = rr[-1]
        return arr

    def _assign_single_chain_fg_cum_extend_type(
        self, molecule, a, fgn, fg, tag, arr, not_a_center
    ):
        atom = molecule.Atoms[a]
        arr.append(a)
        atom.function_group_tag = f"{tag}_center"
        atom.function_group_id = fgn
        _extend = [
            [
                an,
                molecule.Atoms[an].element.lower(),
                self._a_symbol[molecule.Atoms[an].element],
            ]
            for an in atom.connectivity
            if an not in arr + not_a_center and molecule.Atoms[an].element not in ["C"]
        ]
        _extend = sorted(_extend, key=lambda x: x[2])
        fg = f"{fg}{''.join([rr[1] for rr in _extend])}"
        atom.function_group = fg
        _tmp = [an for an in atom.connectivity if an not in arr + not_a_center]
        for an in _tmp:
            arr.append(an)
            molecule.Atoms[an].function_group = fg
            molecule.Atoms[an].function_group_id = fgn
            molecule.Atoms[an].function_group_tag = (
                f"{molecule.Atoms[an].element}_aplpha"
            )
        return arr

    def _assign_single_chain_fg_double_extend_type(
        self, molecule, a, fgn, fg, tag, arr, couple_pair, not_a_center, a_center_dict
    ):

        def _get_connects(atom):
            _extend_ = [
                [-1, "", "z"],
                [-1, "", "z"],
                [-1, "", "z"],
                [-1, "", "z"],
                [-1, "", "z"],
                [-1, "", "z"],
            ]
            for ii, an in enumerate(atom.connectivity):
                if an not in arr + not_a_center and an in a_center_dict:
                    _extend_[ii] = [
                        an,
                        molecule.Atoms[an].element.lower(),
                        self._a_symbol[molecule.Atoms[an].element],
                    ]
            _extend_ = sorted(_extend_, key=lambda x: x[2])
            return [
                [rr[0] for rr in _extend_[:2] if rr[0] != -1],
                f"{_extend_[0][1]}{_extend_[1][1]}",
                f"{_extend_[0][2]}{_extend_[1][2]}",
            ]

        atom = molecule.Atoms[a]
        arr.append(a)
        atom.function_group_tag = f"{tag}_center"
        atom.function_group_id = fgn
        catom = molecule.Atoms[couple_pair[a]]
        arr.append(couple_pair[a])
        catom.function_group_tag = f"{catom.element.lower()}_center"
        catom.function_group_id = fgn
        _extend_0 = _get_connects(atom)
        _extend_1 = _get_connects(catom)
        _extend = sorted([_extend_0, _extend_1], key=lambda x: x[2])
        _fg_order = f"{self._fg_order.index(fg)}:{_extend[0][2]}{_extend[1][2]}"
        if len(_extend[1][0]) == 0:
            _extend = [_extend[1], _extend[0]]
        fg = f"{_extend[0][1]}{fg}{_extend[1][1]}"
        atom.function_group = fg
        catom.function_group = fg
        atom._double_fg_order = _fg_order
        catom._double_fg_order = _fg_order
        for rr in _extend:
            for rrr in rr[0]:
                an = rrr
                arr.append(an)
                molecule.Atoms[an].function_group = fg
                molecule.Atoms[an].function_group_id = fgn
                molecule.Atoms[an].function_group_tag = (
                    f"{molecule.Atoms[an].element}_alpha"
                )
                molecule.Atoms[an]._double_fg_order = _fg_order
        return arr

    def _assign_sing_chain_fg(
        self, molecule, a, fgn, arr, couple_pair, not_a_center, a_center_dict
    ):
        fg, tag = molecule.Atoms[a].atom_type_name.split(":")
        fg = fg.strip("_")
        if fg in self._fgs["bsi"] + self._fgs["acid"] + self._fgs["hacid"]:
            molecule.Atoms[a]._protected_fg = fg
        if fg in self._fgs["type_1"]:
            arr = self._assign_sing_chain_fg_acid_type(
                molecule, a, fgn, fg, tag, arr, not_a_center
            )
        elif fg in self._fgs["type_4"]:
            arr = self._assign_single_chain_fg_b_extend_type(
                molecule, a, fgn, fg, tag, arr, couple_pair, not_a_center, a_center_dict
            )
        elif fg in self._fgs["type_2"]:
            arr = self._assign_single_chain_fg_cum_extend_type(
                molecule, a, fgn, fg, tag, arr, not_a_center
            )
        elif fg in self._fgs["type_5"] + self._fgs["type_6"]:
            arr = self._assign_single_chain_fg_double_extend_type(
                molecule, a, fgn, fg, tag, arr, couple_pair, not_a_center, a_center_dict
            )
        elif fg in self._fgs["type_7"]:
            arr = self._assign_sing_chain_fg_a_type(molecule, a, fgn, fg, tag, arr)
        else:
            arr = self._assign_sing_chain_fg_non_extend_type(
                molecule, a, fgn, fg, tag, arr, not_a_center
            )
        return arr

    def _combind_a(self, molecule):

        def _get_fg(ans, molecule):
            name1 = [molecule.Atoms[an].function_group for an in ans]
            name2 = [molecule.Atoms[an].function_group for an in reversed(ans)]
            label1 = "".join([self._a_symbol[ff] for ff in name1])
            label2 = "".join([self._a_symbol[ff] for ff in name2])
            return (
                sorted(
                    [["".join(name1), label1], ["".join(name2), label2]],
                    key=lambda x: x[1],
                )[0],
                min([molecule.Atoms[an].function_group_id for an in ans]),
            )

        a_connect_dict = {
            atom.ID: atom.connectivity
            for atom in molecule.Atoms
            if getattr(atom, "function_group", None) in self._fgs["a_type"]
            and atom.element not in ["H"]
        }
        for frag in find_mole(a_connect_dict):
            connect_dict = {
                ii: [an for an in molecule.Atoms[ii].connectivity if an in frag]
                for ii in frag
            }
            terminal_atom_arr = [
                ii for ii, conn in connect_dict.items() if len(conn) <= 1
            ]
            all_chain = search_chain(terminal_atom_arr, connect_dict)
            for chain in all_chain:
                _this_fg, _this_id = _get_fg(chain, molecule)
                for an in chain:
                    molecule.Atoms[an].function_group = _this_fg[0]
                    molecule.Atoms[an]._combine_a_fg_order = _this_fg[1]
                    molecule.Atoms[an].function_group_id = _this_id

    def _check_H(self, molecule):
        for atom in molecule.Atoms:
            if atom.element in SCAOG and atom.connect:
                atom_t = molecule.Atoms[atom.connect[0]]
                atom.function_group = atom_t.function_group
                atom.function_group_id = atom_t.function_group_id
                atom.function_group_tag = f"{atom.elem}_alpha"

    def _check_other(self, molecule, fgn):
        this_connectivities = {
            atom.ID: atom.connectivity
            for atom in molecule.Atoms
            if not hasattr(atom, "function_group")
            and (
                atom.element not in SCAOG
                or all(molecule.Atoms[a].element in SCAOG for a in atom.connect)
            )
        }
        for frag in find_mole(this_connectivities):
            fgn += 1
            for an in frag:
                atom = molecule.Atoms[an]
                atom.function_group = "other"
                atom.function_group_id = fgn
                atom.function_group_tag = f"{atom.elem}_{len(atom.connectivity)}"

    def _molecule_chain_fg_order(self, molecule, res_atoms):
        orders = []
        for attribute in ("_a_b_a_fg_order", "_combine_a_fg_order", "_double_fg_order"):
            by_group = {}
            for atom in molecule.Atoms:
                if hasattr(atom, attribute):
                    group, value = atom.function_group, getattr(atom, attribute)
                    # A terminal atom may precede the metadata carrier. Inspect
                    # every carrier and resolve repeated labels deterministically.
                    if group not in by_group or value < by_group[group]:
                        by_group[group] = value
            orders.append(sorted(by_group, key=lambda group: (by_group[group], group)))
        return tuple(orders)

    def _anhy_fg(self, molecule, acid_fg, _a_b_a_fg):

        def _judge_fg(atom, molecule, _this_fg, center=False, alpha=False):
            if atom.function_group not in _this_fg:
                return False
            if not alpha:
                if atom.function_group_tag.find("center") != -1:
                    return True
                elif center:
                    return False
            __nn = [
                atom.bond_type[ii]
                for ii, an in enumerate(atom.connectivity)
                if molecule.Atoms[an].function_group_tag.find("center") != -1
            ]
            if len(__nn) > 0 and "2" not in __nn:
                return True
            return False

        def _set_anhy_fg(molecule, fgns, fgs, fgs_reverse, _flag=False):
            if not _flag:
                new_fg = "_".join(fgs)
                new_fg_reverse = "_".join(fgs_reverse)
            else:
                new_fg = "".join(fgs)
                new_fg_reverse = "".join(fgs_reverse)
            fgn = min(fgns)
            for atom in molecule.Atoms:
                if atom.function_group_id in fgns:
                    atom._function_group_direction = (
                        "left" if atom.function_group_id == fgns[0] else "right"
                    )
                    atom.function_group_id = fgn
                    atom.function_group = new_fg
                    atom._function_group_reverse = new_fg_reverse
            return new_fg

        def _get_side_fg(atom, reverse=False, _acid_flag=False):

            def _get_tag(atom):
                _tmp_ = atom.atom_type_name.split(":")
                if len(_tmp_) == 2:
                    _tag = f"{_tmp_[0]}"
                else:
                    _tag = f"{atom.element}d"
                return _tag

            if _acid_flag:
                _tag = _get_tag(atom)
            else:
                _tag = ""
            _direction = getattr(atom, "_function_group_direction", None)
            if _direction is None:
                atom_function_group = atom.function_group
                atom_function_group_reverse = atom.function_group
            elif _direction == "left":
                atom_function_group = (
                    atom._function_group_reverse if not reverse else atom.function_group
                )
                atom_function_group_reverse = (
                    atom._function_group_reverse if reverse else atom.function_group
                )
            elif _direction == "right":
                atom_function_group = (
                    atom.function_group if not reverse else atom._function_group_reverse
                )
                atom_function_group_reverse = (
                    atom.function_group if reverse else atom._function_group_reverse
                )
            if not reverse:
                atom_function_group_reverse = (
                    atom_function_group_reverse
                    if not _acid_flag
                    else f"{_tag}_{atom_function_group_reverse}"
                )
                atom_function_group = (
                    atom_function_group
                    if not _acid_flag
                    else f"{atom_function_group}_{_tag}"
                )
            else:
                atom_function_group_reverse = (
                    atom_function_group_reverse
                    if not _acid_flag
                    else f"{atom_function_group_reverse}_{_tag}"
                )
                atom_function_group = (
                    atom_function_group
                    if not _acid_flag
                    else f"{_tag}_{atom_function_group}"
                )
            return (atom_function_group, atom_function_group_reverse)

        _this_fg = acid_fg + _a_b_a_fg
        _candinate_anhy_connect_atoms = [
            [atom, _this_fg.index(atom.function_group)]
            for atom in molecule.Atoms
            if _judge_fg(atom, molecule, _this_fg)
        ]
        _candinate_anhy_connect_atoms = sorted(
            _candinate_anhy_connect_atoms, key=lambda x: x[1]
        )
        _candinate_anhy_connect_atoms_arr = [
            rr[0].ID for rr in _candinate_anhy_connect_atoms
        ]
        _candinate_anhy_target_alpha_atoms_arr = [
            rr[0].ID
            for rr in _candinate_anhy_connect_atoms
            if _judge_fg(rr[0], molecule, _this_fg, center=True)
        ]
        _candinate_anhy_target_alpha_atoms_arr = list(
            set(_candinate_anhy_connect_atoms_arr).difference(
                set(_candinate_anhy_target_alpha_atoms_arr)
            )
        )
        _candinate_anhy_target_acid_alpha_atoms_arr = [
            an
            for an in _candinate_anhy_target_alpha_atoms_arr
            if molecule.Atoms[an].function_group in acid_fg
        ]
        _ane_atoms_ = [
            atom
            for atom in molecule.Atoms
            if atom.function_group in self._fgs["type_8"] and atom.element not in SCAOG
        ]
        _center_atoms_ = [
            [atom.ID, _this_fg.index(atom.function_group)]
            for atom in molecule.Atoms
            if _judge_fg(atom, molecule, _this_fg, center=True)
        ]
        _center_atoms_ = sorted(_center_atoms_, key=lambda x: x[1])
        _center_atoms_arr = [rr[0] for rr in _center_atoms_]
        _anhy_fg = []
        for an in _candinate_anhy_connect_atoms_arr:
            atom = molecule.Atoms[an]
            _conn_fgn = [
                molecule.Atoms[an]
                for an in atom.connectivity
                if an in _candinate_anhy_connect_atoms_arr
                and molecule.Atoms[an].function_group_id != atom.function_group_id
            ]
            if len(_conn_fgn) > 0:
                atom.function_group_tag += "_canhy"
                _flag = False
                _c_acid_flag = False
                _acid_flag = False
                if _conn_fgn[0].ID in _candinate_anhy_target_alpha_atoms_arr:
                    _flag = True
                if _conn_fgn[0].ID in _candinate_anhy_target_acid_alpha_atoms_arr:
                    _c_acid_flag = True
                if atom.ID in _candinate_anhy_target_acid_alpha_atoms_arr:
                    _acid_flag = True
                if not _flag:
                    catom_function_group = _conn_fgn[0].function_group
                    catom_function_group_reverse = getattr(
                        _conn_fgn[0],
                        "_function_group_reverse",
                        _conn_fgn[0].function_group,
                    )
                else:
                    catom_function_group, catom_function_group_reverse = _get_side_fg(
                        _conn_fgn[0], reverse=True, _acid_flag=_c_acid_flag
                    )
                atom_function_group, atom_function_group_reverse = _get_side_fg(
                    atom, _acid_flag=_acid_flag
                )
                fgns = [atom.function_group_id, _conn_fgn[0].function_group_id]
                fgs = [atom_function_group, catom_function_group]
                fgs_reverse = [
                    catom_function_group_reverse,
                    atom_function_group_reverse,
                ]
                _anhy_fg.append(
                    _set_anhy_fg(molecule, fgns, fgs, fgs_reverse, _flag=_flag)
                )
        tmp = [[fg, len(fg.split("_"))] for fg in _anhy_fg]
        tmp = sorted(tmp, key=lambda x: x[1], reverse=True)
        _ane_anhy_fg = []
        for atom in _ane_atoms_:
            _conn_fgn = [
                molecule.Atoms[an]
                for an in atom.connectivity
                if an in _center_atoms_arr
            ]
            if len(_conn_fgn) >= 2:
                _conn_fgn = [
                    molecule.Atoms[rr[0]]
                    for rr in sorted(
                        [
                            [catom.ID, _center_atoms_arr.index(catom.ID)]
                            for catom in _conn_fgn
                        ],
                        key=lambda x: x[1],
                    )
                ]
                atom.function_group_tag += "_cacid"
                atom0, atom1 = _conn_fgn[:2]
                atom0_function_group, atom0_function_group_reverse = _get_side_fg(atom0)
                atom1_function_group, atom1_function_group_reverse = _get_side_fg(
                    atom1, reverse=True
                )
                fgns = [
                    atom0.function_group_id,
                    atom.function_group_id,
                    atom1.function_group_id,
                ]
                fgs = [atom0_function_group, atom.function_group, atom1_function_group]
                fgs_reverse = [
                    atom0_function_group_reverse,
                    atom.function_group,
                    atom1_function_group_reverse,
                ]
                _ane_anhy_fg.append(_set_anhy_fg(molecule, fgns, fgs, fgs_reverse))
        tmp_ane = [[fg, len(fg.split("_"))] for fg in _ane_anhy_fg]
        tmp_ane = sorted(tmp_ane, key=lambda x: x[1], reverse=True)
        return [fg[0] for fg in tmp] + [fg[0] for fg in tmp_ane]

    def _assign_chain_fg(self, molecule, fgn, res_atoms, arr):
        fgn += 1
        fg_center, not_a_center, couple_pair, a_center_dict = self._get_chain_center(
            molecule, res_atoms
        )
        arr.extend([vv for vv in couple_pair.values()])
        for ii, rr in enumerate(fg_center):
            an = rr[0]
            if an not in arr:
                fgn += 1
                arr = self._assign_sing_chain_fg(
                    molecule, an, fgn, arr, couple_pair, not_a_center, a_center_dict
                )
        self._combind_a(molecule)
        self._check_other(molecule, len(fg_center) + fgn)
        self._check_H(molecule)
        _a_b_a_fg, _combine_a_fg, _double_fg = self._molecule_chain_fg_order(
            molecule, res_atoms
        )
        _anhy_fg = []
        _anhy_fg = self._anhy_fg(
            molecule, self._fgs["type_1"] + self._fgs["bsi"], _a_b_a_fg + self._fgs["b"]
        )
        return (_a_b_a_fg, _combine_a_fg, _double_fg, _anhy_fg)

    def _assign_fused_ring_fg(self, molecule, fgn):
        arr = []
        small_ring_atoms = []
        for __, ring in molecule.rings.items():
            if len(ring[:-1]) < 10:
                small_ring_atoms.extend(ring[:-1])
        small_ring_atoms = set(small_ring_atoms)
        _this_connectivity = {
            ai: [
                an
                for ii, an in enumerate(molecule.Atoms[ai].connectivity)
                if an in small_ring_atoms
                and molecule.Atoms[ai].connectivity_type[ii] not in BREAK_BOND
            ]
            for ai in small_ring_atoms
        }
        _tmp = find_mole(_this_connectivity)
        for rr in _tmp:
            fgn += 1
            for an in rr:
                molecule.Atoms[an].function_group = "R"
                molecule.Atoms[an].function_group_id = fgn
                molecule.Atoms[an].function_group_tag = (
                    f"r_{molecule.Atoms[an].element.lower()}_{molecule.Atoms[an].local}"
                )
                arr.append(an)
        return (arr, {"R": -1}, fgn)

    def _assign_nonfused_ring_fg(self, molecule, fgn):
        arr = []
        fgs = {}
        rings = deepcopy(molecule._rings)
        special_ring_atoms = []
        total_ring_atoms = []
        for __, ring in rings.items():
            rn = len(ring[:-1])
            if rn in [3, 4, 5, 6, 7, 8, 9]:
                special_ring_atoms.extend(ring[:-1])
                fg = f"{len(ring[:-1])}R"
                elems = [molecule.Atoms[ii].elem for ii in ring[:-1]]
                for ee in ["O", "N", "S", "P"]:
                    if ee in elems:
                        fg += ee
                if ring[-1] in ["ar1", "ar2", "ar3", "ar4", "ar5"]:
                    fg += f"_{ring[-1]}"
                ring.append(fg)
                ring.append(fgn)
                if rn > 6:
                    _mfg = "MR"
                else:
                    _mfg = fg
                if _mfg not in fgs:
                    fgs[_mfg] = []
                fgn += 1
                fgs[_mfg].append(fgn)
            else:
                fgs["LR"] = [-1]
            total_ring_atoms.extend(ring[:-3])
        special_ring_atoms = list(set(special_ring_atoms))
        for ii in special_ring_atoms:
            atom = molecule.Atoms[ii]
            _order = ["ar1", "ar2", "ar3", "ar4", "ar5"]
            tmp = [[], [], [], [], [], []]
            for nn, ar in enumerate(atom.has_ring_property):
                if ar in _order:
                    index_jj = _order.index(ar)
                    tmp[index_jj].append([atom.has_ring_size[nn], atom.has_ring[nn]])
                else:
                    tmp[5].append([atom.has_ring_size[nn], atom.has_ring[nn]])
            for rr in tmp:
                if len(rr) != 0:
                    atom.function_group = rings[sorted(rr, key=lambda x: x[0])[0][1]][
                        -2
                    ]
                    atom.function_group_id = rings[
                        sorted(rr, key=lambda x: x[0])[0][1]
                    ][-1]
                    break
            tag = f"{atom.elem.lower()}_{atom.local}"
            if len(set(atom.has_ring_property).intersection({"ar1", "ar2"})) > 0:
                bts = [
                    atom.bond_type[jj]
                    for jj, an in enumerate(atom.connect)
                    if an in total_ring_atoms
                ]
                if set(bts) == set(["1"]):
                    tag += "_pc"
                else:
                    tag += "_pi"
            atom.function_group_tag = tag
            arr.append(ii)
        return (arr, fgs, fgn)

    def _assign_ring_fg(self, molecule, fgn=-1):
        arr, fgs, fgn = self._assign_fused_ring_fg(molecule, fgn)
        special_ring_atoms = deepcopy(arr)
        for ii in special_ring_atoms:
            atom = molecule.Atoms[ii]
            for nn, bt in enumerate(atom.bond_type):
                an1 = atom.connect[nn]
                if an1 not in special_ring_atoms:
                    atom1 = molecule.Atoms[an1]
                    if bt == "2":
                        atom1.function_group = atom.function_group
                        atom1.function_group_id = atom.function_group_id
                        atom1.function_group_tag = (
                            f"{atom1.elem.lower()}_{atom1.local}_exdb"
                        )
                        arr.append(an1)
                        for jj, an2 in enumerate(atom1.connectivity):
                            atom2 = molecule.Atoms[an2]
                            if atom1.bond_type[jj] != "1" and atom2.local in LOCAL_EN:
                                atom2.function_group = atom.function_group
                                atom2.function_group_id = atom.function_group_id
                                atom2.function_group_tag = (
                                    f"{atom1.elem.lower()}_beta_exdb"
                                )
                                arr.append(an2)
                        atom.function_group_tag += "_cexdb"
                    elif atom.formal_charge != 0 and atom1.formal_charge != 0:
                        atom1.function_group = atom.function_group
                        atom1.function_group_id = atom.function_group_id
                        atom1.function_group_tag = (
                            f"{atom1.elem.lower()}_{atom1.local}_ion"
                        )
                        arr.append(an1)
                        atom.function_group_tag += "_cion"
                    elif atom1.elem in SCAOG:
                        atom1.function_group = atom.function_group
                        atom1.function_group_id = atom.function_group_id
                        atom1.function_group_tag = f"{atom1.elem.lower()}_{atom1.local}"
                        arr.append(an1)
                    else:
                        if atom1.element == "C":
                            _H_atoms = [
                                ann
                                for ann in atom1.connectivity
                                if molecule.Atoms[ann].element in SCAOG
                            ]
                        else:
                            _H_atoms = [
                                ann
                                for ann in atom1.connectivity
                                if molecule.Atoms[ann].element in ["H"]
                            ]
                        if (
                            len(_H_atoms) - atom1.formal_charge
                            == Element.get(atom1.element).valents[0] - 1
                        ):
                            atom1.function_group = atom.function_group
                            atom1.function_group_id = atom.function_group_id
                            atom1.function_group_tag = (
                                f"{atom1.elem.lower()}_{atom1.local}"
                            )
                            arr.append(an1)
                            for ann in _H_atoms:
                                molecule.Atoms[ann].function_group = atom.function_group
                                molecule.Atoms[ann].function_group_id = (
                                    atom.function_group_id
                                )
                                molecule.Atoms[ann].function_group_tag = (
                                    atom1.function_group_tag + "_alpha"
                                )
        _tmp = []
        for kk, vv in fgs.items():
            if kk == "LR" or kk == "R":
                _tmp.append([kk, vv, -2])
            elif kk == "MR":
                _tmp.append([kk, vv, -1])
            else:
                _tmp.append([kk, vv, self._fgs["rings"].index(kk)])
        _tmp = sorted(_tmp, key=lambda x: x[2])
        return (arr, [[rr[0], rr[1]] for rr in _tmp], fgn)

    def resort_fg_id(self, molecule):
        fgs = {}
        for atom in molecule.Atoms:
            if atom.function_group_id not in fgs:
                fgs[atom.function_group_id] = []
            fgs[atom.function_group_id].append(atom.ID)
        nn = -1
        tmp = []
        for atom in molecule.Atoms:
            if atom.ID not in tmp:
                nn += 1
                tmp.extend(fgs[atom.function_group_id])
                for an in fgs[atom.function_group_id]:
                    molecule.Atoms[an].function_group_id = nn

    def _molecule_fg_order(self, molecule):
        this_order = []
        for nn in range(0, 9):
            if nn == 1:
                this_order.extend(molecule._anhy_fgs)
            elif nn == 4:
                this_order.extend(molecule._a_b_a_fgs)
            elif nn == 5:
                this_order.extend(molecule._double_fgs)
            elif nn == 7:
                this_order.extend(molecule._combine_a_fgs)
            elif nn == 8:
                this_order.extend(["other"])
            this_order.extend(self._fgs[f"type_{nn}"])
        _tmp = {}
        for atom in molecule.Atoms:
            if atom.function_group not in self._fgs["rings"]:
                if atom.function_group not in _tmp:
                    _tmp[atom.function_group] = [
                        atom.function_group,
                        this_order.index(atom.function_group),
                        [],
                    ]
                _tmp[atom.function_group][2].append(atom.function_group_id)
        _total_fgs = sorted([vv for vv in _tmp.values()], key=lambda x: x[1])
        return [[rr[0], list(set(rr[2]))] for rr in _total_fgs]

    def _combine_connect_ane(self, molecule):

        def _get_connect_order(molecule):
            _connect_order = []
            _tmp_fg = [
                molecule._anhy_fgs
                + self._fgs["inorganic"]
                + self._fgs["acid"]
                + self._fgs["hacid"]
                + self._fgs["ion1"]
                + self._fgs["bsi"],
                self._fgs["acid"]
                + self._fgs["hacid"]
                + molecule._anhy_fgs
                + molecule._a_b_a_fgs,
                self._fgs["a1"] + self._fgs["ion2"] + self._fgs["a"],
                self._fgs["ion1"] + self._fgs["bsi"] + molecule._double_fgs,
                molecule._a_b_a_fgs,
                molecule._combine_a_fgs + self._fgs["cum"],
                molecule._anhy_fgs + molecule._a_b_a_fgs + self._fgs["b"],
                molecule._double_fgs,
                self._fgs["tribond"]
                + self._fgs["az"]
                + self._fgs["yne"]
                + self._fgs["ene"],
            ]
            _fg_order = [f"{r}:center" for r in _tmp_fg[0]]
            _fg_order.extend([f"{r}:side" for r in _tmp_fg[1]])
            _fg_order.extend([f"{r}:center" for r in _tmp_fg[2]])
            _fg_order.extend([f"{r}:side" for r in _tmp_fg[3]])
            _fg_order.extend([f"{r}:2-side" for r in _tmp_fg[4]])
            _fg_order.extend([f"{r}:side" for r in _tmp_fg[5]])
            _fg_order.extend([f"{r}:Ccenter" for r in _tmp_fg[6]])
            _fg_order.extend([f"{r}:center" for r in _tmp_fg[7]])
            for atom in molecule.Atoms:
                fg = atom.function_group
                if (
                    atom.element not in SCAOG
                    and fg != "ane"
                    and (
                        len(
                            [
                                an
                                for an in atom.connectivity
                                if molecule.Atoms[an].function_group
                                not in ["ane"] + self._fgs["rings"]
                            ]
                        )
                        > 0
                    )
                ):
                    if "center" in atom.function_group_tag:
                        if fg in _tmp_fg[6]:
                            tag = f"{fg}:Ccenter"
                        else:
                            tag = f"{fg}:center"
                    elif fg in _tmp_fg[4] and "beta" in atom.function_group_tag:
                        tag = f"{fg}:2-side"
                    else:
                        tag = f"{fg}:side"
                    if tag in _fg_order:
                        _connect_order.append([atom.ID, _fg_order.index(tag)])
            _connect_order = sorted(_connect_order, key=lambda x: x[1])
            return [rr[0] for rr in _connect_order]

        _connect_order = _get_connect_order(molecule)
        res_atom = []
        for atom in molecule.Atoms:
            if atom.element in SCAOG:
                continue
            if atom.function_group == "ane" and atom.element == "C":
                _canidate = []
                _tmp = [atom.ID]
                for an in atom.connectivity:
                    if molecule.Atoms[an].function_group_id == atom.function_group_id:
                        _tmp.append(an)
                    elif an in _connect_order:
                        _canidate.append([an, _connect_order.index(an)])
                if len(_canidate) > 0:
                    connect_an = sorted(_canidate, key=lambda x: x[1])[0][0]
                    _this_function_group = molecule.Atoms[connect_an].function_group
                    _this_function_group_id = molecule.Atoms[
                        connect_an
                    ].function_group_id
                    _this_function_group_tag = "ane_extend"
                    for an in _tmp:
                        molecule.Atoms[an].function_group = _this_function_group
                        molecule.Atoms[an].function_group_id = _this_function_group_id
                        molecule.Atoms[an].function_group_tag = (
                            f"{_this_function_group_tag}_{molecule.Atoms[an].element.lower()}"
                        )
                else:
                    res_atom.append(atom.ID)
        a_atom = []
        cum_atom = []
        for atom in molecule.Atoms:
            if atom.function_group in self._fgs["a_type"] and atom.element not in SCAOG:
                _ss = [
                    an
                    for an in atom.connectivity
                    if molecule.Atoms[an].function_group_id == atom.function_group_id
                    and molecule.Atoms[an].element not in SCAOG
                ]
                if len(_ss) == 0:
                    a_atom.append(atom.ID)
            elif (
                atom.function_group in self._fgs["type_2"] and atom.element not in SCAOG
            ):
                _ss = [
                    an
                    for an in atom.connectivity
                    if molecule.Atoms[an].function_group_id == atom.function_group_id
                    and molecule.Atoms[an].element not in SCAOG
                ]
                if len(_ss) == 0:
                    cum_atom.append(atom.ID)
        return (res_atom, a_atom, cum_atom)

    def _linear_chain(self, rr):
        arr = rr
        rrr = []
        if len(rr) % 3 == 2:
            rrr.append([rr[0], rr[1]])
            arr = rr[2:]
        elif len(rr) % 3 == 1:
            rrr.append([rr[0], rr[1]])
            rrr.append([rr[-2], rr[-1]])
            arr = rr[2:-2]
        if len(arr) % 3 == 0:
            for i in range(0, len(arr), 3):
                rrr.append([arr[i], arr[i + 1], arr[i + 2]])
        return rrr

    def _ane_chain(self, molecule, res_atom):
        _connectivity = {
            ai: [an for an in molecule.Atoms[ai].connectivity if an in res_atom]
            for ai in res_atom
        }
        _tmp = find_mole(_connectivity)
        _tmp_tmp = []
        _branch_tmp = []
        for rr in _tmp:
            if len(rr) < 4:
                _tmp_tmp.append(rr)
            else:
                _branch_tmp.append(rr)
        if len(_branch_tmp) > 0:
            for rr in _branch_tmp:
                _connect_dict = {
                    rn: sorted(
                        [an for an in molecule.Atoms[rn].connect if an in rr],
                        key=lambda a: molecule._stable_atom_ranks[a],
                    )
                    for rn in sorted(rr, key=lambda a: molecule._stable_atom_ranks[a])
                }
                _terminal_atom = [
                    rn for rn in _connect_dict if len(_connect_dict[rn]) == 1
                ]
                branchs = search_chain(_terminal_atom, _connect_dict)
                if len(branchs[0]) < 4:
                    _tmp_tmp.append(branchs[0])
                else:
                    _tmp_tmp.extend(self._linear_chain(branchs[0]))
                for branch in branchs[1:]:
                    if len(branch[1:]) < 4:
                        _tmp_tmp.append(branch[1:])
                    else:
                        _tmp_tmp.extend(self._linear_chain(branch[1:]))
        single_atom = []
        for rr in _tmp_tmp:
            if len(rr) == 1:
                single_atom.append(rr[0])
            else:
                fgn = min([molecule.Atoms[r].function_group_id for r in rr])
                for r in rr:
                    r_fgn = molecule.Atoms[r].function_group_id
                    for an in molecule.Atoms[r].connectivity:
                        if molecule.Atoms[an].function_group_id == r_fgn:
                            molecule.Atoms[an].function_group_id = fgn
                    molecule.Atoms[r].function_group_id = fgn
        total_tmp = []
        for rr in _tmp_tmp:
            total_tmp.extend(rr)
        single_atom.extend(list(set(res_atom).difference(set(total_tmp))))
        return single_atom

    def _sing_atom(self, molecule, single_atom, break_order):
        for an in sorted(single_atom, key=lambda a: molecule._stable_atom_ranks[a]):
            _t = [an]
            _tmp = []
            for ann in molecule.Atoms[an].connectivity:
                if (
                    molecule.Atoms[ann].function_group_id
                    == molecule.Atoms[an].function_group_id
                ):
                    _t.append(ann)
                elif molecule.Atoms[ann].function_group in break_order:
                    if molecule.Atoms[ann].function_group_tag.find("ane_extend") != -1:
                        _tmp.append(
                            [
                                ann,
                                len(break_order)
                                + break_order.index(molecule.Atoms[ann].function_group),
                            ]
                        )
                    else:
                        _tmp.append(
                            [ann, break_order.index(molecule.Atoms[ann].function_group)]
                        )
            if len(_tmp) > 0:
                ann = min(
                    _tmp, key=lambda x: (x[1], molecule._stable_atom_ranks[x[0]])
                )[0]
                fgn = molecule.Atoms[ann].function_group_id
                fg = molecule.Atoms[ann].function_group
                tag = f"{molecule.Atoms[an].function_group}_extend"
                for a in _t:
                    molecule.Atoms[a].function_group_id = fgn
                    molecule.Atoms[a].function_group = fg
                    molecule.Atoms[a].function_group_tag = tag

    def _combine_singe_atoms(self, molecule, res_atom, a_atom, cum_atom, break_order):
        single_cum_atom = self._ane_chain(molecule, cum_atom)
        single_ane_atom = self._ane_chain(molecule, res_atom)
        single_atom = a_atom + single_cum_atom + single_ane_atom
        if len(single_atom) > 0:
            self._sing_atom(molecule, single_atom, break_order)

    def _molecule_split(self, molecule):
        for atom in molecule.Atoms:
            atom.break_bond = [None for _ in atom.connectivity]
            for ii, an in enumerate(atom.connectivity):
                atom.break_bond[ii] = (
                    False
                    if atom.function_group_id == molecule.Atoms[an].function_group_id
                    else True
                )
        _connectivity = {
            atom.ID: [
                an for ii, an in enumerate(atom.connectivity) if not atom.break_bond[ii]
            ]
            for atom in molecule.Atoms
        }
        return find_mole(_connectivity)

    def _save_double_bond(self, molecule):
        for bond in molecule.Bonds:
            atom0 = molecule.Atoms[bond.a1]
            atom1 = molecule.Atoms[bond.a2]
            if atom0.function_group_id != atom1.function_group_id:
                flag = False
                if bond.get_type(molecule) != "1":
                    flag = True
                if (
                    hasattr(atom0, "_protected_fg") or hasattr(atom1, "_protected_fg")
                ) and (
                    atom0.local not in LOCAL_CYCLE and atom1.local not in LOCAL_CYCLE
                ):
                    flag = True
                if flag:
                    fgn = min([atom0.function_group_id, atom1.function_group_id])
                    ans = [
                        atom
                        for atom in molecule.Atoms
                        if atom.function_group_id
                        in [atom0.function_group_id, atom1.function_group_id]
                    ]
                    for atom in ans:
                        atom.function_group_id = fgn

    def run(self, molecule):
        if "fragmentation" in molecule.steps:
            self._clean(molecule)
        tmp, _ring_fgs, rfgn = self._assign_ring_fg(molecule)
        res_atoms = list(set((atom.ID for atom in molecule.Atoms)).difference(set(tmp)))
        self._assign_atom_type(molecule, atom_arr=res_atoms)
        (
            molecule._a_b_a_fgs,
            molecule._combine_a_fgs,
            molecule._double_fgs,
            molecule._anhy_fgs,
        ) = self._assign_chain_fg(molecule, rfgn, res_atoms, tmp)
        res_ane_atom, a_atom, cum_atom = self._combine_connect_ane(molecule)
        total_fgs = self._molecule_fg_order(molecule)
        self._combine_singe_atoms(
            molecule, res_ane_atom, a_atom, cum_atom, [rr[0] for rr in total_fgs]
        )
        self.resort_fg_id(molecule)
        self._save_double_bond(molecule)
        molecule.frags = self._molecule_split(molecule)
        molecule.steps.append("fragmentation")
        return molecule
