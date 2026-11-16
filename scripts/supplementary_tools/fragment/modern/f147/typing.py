"""JSON atom-pattern matching for f147 fragmentation only.

Force-field assignment, decorated force-field types and legacy TXT parsing do
not participate in fragmentation. Keep only the pattern tree and type aliases.
"""

import functools
import itertools
import json
import os
from copy import deepcopy


class TypingDefine:
    def __init__(self, atversion, check_at_types=True):
        with open(atversion) as stream:
            raw = json.load(stream)
        self.at_types = {
            name: {key: rule[key] for key in ("children", "At_Attrs", "At_Connect")}
            for name, rule in raw.items()
        }
        aliases = dict(raw["root"].get("same", {}))
        aliases.update(
            {
                name: rule["same"]
                for name, rule in raw.items()
                if name != "root" and rule.get("same")
            }
        )
        self.at_types["root"]["same"] = aliases
        if check_at_types:
            self._check_at_types()

    def _check_at_types(self):
        for name, rule in self.at_types.items():
            if set(rule["children"]) - self.at_types.keys():
                raise ValueError(f"Undefined child type in {name}")
            for attrs in rule["At_Attrs"].values():
                for key, values in attrs.items():
                    expected = (
                        int
                        if key in ("coor", "ring")
                        else (int, float) if key == "charge" else None
                    )
                    if expected and any(
                        value != "?" and not isinstance(value, expected)
                        for value in values
                    ):
                        raise ValueError(f"Invalid {key} values in {name}")


@functools.lru_cache(maxsize=8)
def _compiled_rules(path, check_at_types, file_signature):
    """Private parsed template; callers copy it before molecular matching.

    The file signature is part of the cache key so edits invalidate the entry.
    No molecule, matching state or typing engine is cached.
    """
    return TypingDefine(path, check_at_types=check_at_types).at_types


class TypingEngine:
    def __init__(self, atversion, extend_depth=8, check_at_types=True):
        self._extend_depth = extend_depth
        path = os.path.realpath(atversion)
        stat = os.stat(path)
        signature = (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size)
        self.at_types = deepcopy(_compiled_rules(path, check_at_types, signature))

    def assign_mole_at(self, molecule, atoms_arr=None):
        if "atom type" in molecule.steps:
            return molecule
        if atoms_arr is None:
            atoms_arr = range(len(molecule.Atoms))
        for a in atoms_arr:
            self.search_at_tree(a, molecule.Atoms)
            molecule.Atoms[a].ats_tree = self.ats_tree
            atom_type = self.ats_tree[-1]
            molecule.Atoms[a].atom_type_name = self.at_types["root"]["same"].get(
                atom_type, atom_type
            )
        return molecule

    def search_at_tree(self, a, Atoms):
        self.ats_tree, self.find_flag = ([], False)
        for at in self.at_types["root"]["children"]:
            self.find_at(a, at, Atoms)
            if len(self.ats_tree) != 0:
                break

    def find_at(self, a, at, Atoms):
        if not self.find_flag and self.match_at_a(at, a, Atoms):
            self.ats_tree.append(at)
            for at0 in self.at_types[at]["children"]:
                self.find_at(a, at0, Atoms)
            self.find_flag = True

    def match_at_a(self, at, a, Atoms):
        this_connect = self.at_types[at]["At_Connect"]
        this_attrs = self.at_types[at]["At_Attrs"]
        if not self.match_at_attrs(Atoms[a], this_attrs["0"]):
            return False
        match_atom = [[] for i in range(len(this_attrs))]
        match_atom[0].append(a)

        def get_candidate_atom(a, n, Atoms):
            arr = [a]
            tmp = set()
            for i in range(0, n):
                res = set(arr) - tmp
                tmp = set(arr)
                for aa in res:
                    arr += getattr(Atoms[aa], "connect")
            return set(arr)

        candidate_atom = get_candidate_atom(a, self._extend_depth, Atoms)
        for i in range(1, len(this_attrs)):
            for aa in candidate_atom:
                if self.match_at_attrs(Atoms[aa], this_attrs[str(i)]):
                    match_atom[i].append(aa)
            if len(match_atom[i]) == 0:
                return False
        for arr in itertools.product(*match_atom):
            if len(set(arr)) != len(this_attrs):
                continue
            for conn in this_connect:
                aa, bb, cc = (arr[conn[0]], arr[conn[1]], conn[2])
                if bb not in Atoms[aa].connect or (
                    cc != "-1"
                    and Atoms[aa].bond_type_aromatic[Atoms[aa].connect.index(bb)] != cc
                ):
                    break
            else:
                return True
        return False

    def match_at_attrs(self, atom, attrs):
        for attr, value in attrs.items():
            if "?" in value:
                continue
            if attr == "coor":
                if len(atom.connect) not in value:
                    return False
            elif attr == "formal_charge":
                if atom.formal_charge not in value:
                    return False
            elif attr == "charge":
                __flag = False
                for cc in value:
                    if abs(atom.primitive_formal_charge - cc) <= 0.01:
                        __flag = True
                        break
                if not __flag:
                    return False
            elif attr == "elem":
                if atom.elem not in value:
                    return False
            elif attr == "ring":
                if len(atom.ring_size) == 0 or min(atom.ring_size) not in value:
                    return False
            elif attr == "arom":
                if all((prop not in value for prop in atom.ring_prop)):
                    return False
            elif attr == "conjugation":
                if not set(atom.bond_type_conjugate).intersection(
                    {"J", "eJ", "fj", "bJ", "j"}
                ):
                    return False
            elif attr == "local":
                if atom.local not in value:
                    return False
            elif attr == "res":
                if atom.residu not in value:
                    return False
            elif attr == "hybrid":
                if atom.hybrid not in value:
                    return False
            elif attr == "btype":
                if len(set(atom.connectivity_type) & set(value)) == 0:
                    return False
        return True


class AtomType:

    def __init__(self, atf, check_at_types=False):
        self.at = TypingEngine(atversion=atf, check_at_types=check_at_types)

    def _atom_type(self, molecules, atoms_arr=None):
        self.at.assign_mole_at(molecules, atoms_arr=atoms_arr)
        aliases = {
            "na+": "Na+",
            "k+": "K+",
            "ca2+": "Ca2+",
            "mg2+": "Mg2+",
            "f-": "F-",
            "cl-": "Cl-",
            "br-": "Br-",
            "i-": "I-",
            "fh": "HF",
            "clh": "HCl",
            "brh": "HBr",
            "ih": "HI",
        }
        for atom in molecules.Atoms:
            prefix, separator, suffix = atom.atom_type_name.partition(":")
            if separator and prefix in aliases:
                atom.atom_type_name = aliases[prefix] + separator + suffix
        if "atom type" not in molecules.steps:
            molecules.steps.append("atom type")
