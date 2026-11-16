"""147 family-generation kernel; product construction belongs to chemistry."""

from copy import deepcopy
import logging
import time
from ..chemistry.chain import search_chain
from ..chemistry.constants import BREAK_BOND, SINGLE_CONNECT_ATOM_OR_GROUP as SCAOG

logger = logging.getLogger(__name__)


class Fragmentation:

    def __init__(self, style="normal"):
        self.style = style

    def run(self, molecule, families):
        """Generate requested families in the established dependency order.

        RSF also marks elementary linker labels consumed by later families and
        exposed in elem_frag. Preserve that side effect even without RSF output.
        CSF must precede those labels, and its visited sets suppress RSF duplicates.
        """
        selected = set(families)
        if not selected:
            return molecule
        self.frags = molecule.frags
        # SSSR uses the prepared parent graph, not elementary descriptors.
        if selected - {"sssr_frag"}:
            self.elementary_frags(molecule)
            if "seco_frag" in selected:
                self.secondary_frags(molecule)
            visited = None
            if selected & {"csf_frag", "rsf_frag"}:
                visited = self.csf_frags(molecule)
            if selected & {
                "elem_frag",
                "rsf_frag",
                "chain_frag",
                "scaffold_frag",
                "tf_frag",
                "sketch_frag",
            }:
                self.rsf_frags(
                    molecule,
                    tmp_frag_index=visited,
                    build_products="rsf_frag" in selected,
                )
            if selected & {"chain_frag", "scaffold_frag", "tf_frag", "sketch_frag"}:
                self.scaffold_tf_frags(
                    molecule,
                    include_scaffold=bool(selected & {"scaffold_frag", "tf_frag"}),
                    include_tf="tf_frag" in selected,
                )
        if "sssr_frag" in selected:
            self.sssr_frags(molecule)
        if "sketch_frag" in selected:
            self.sketch_frags(molecule)
        return molecule

    def elem_frag_label(self, m, frag):
        large_flag_arr = []
        __arom_number = {
            "nonar": 0,
            "ar1": 1,
            "ar2": 2,
            "ar3": 3,
            "ar4": 4,
            "ar5": 5,
            "ar6": 6,
            "ar7": 7,
            "ar8": 8,
        }
        rings = []
        ring_p = []
        atoms = list(
            set(frag).difference(
                set(
                    [
                        aa
                        for aa in frag
                        if m.Atoms[aa].elem in ["H", "F", "Cl", "Br", "I"]
                        and len(m.Atoms[aa].connect) == 1
                    ]
                )
            )
        )
        heavy_atoms = len(atoms)
        ring_atoms = 0
        for aa in frag:
            if len(m.Atoms[aa].inring) > 0:
                large_flag_arr = []
                for rs in m.Atoms[aa].ring_size:
                    if rs >= 10:
                        large_flag_arr.append(True)
                    else:
                        large_flag_arr.append(False)
                if set(large_flag_arr) != {True}:
                    ring_atoms += 1
                    rings.extend(m.Atoms[aa].inring)
                    ring_p.extend(m.Atoms[aa].ring_prop)
        out_ring_atoms = heavy_atoms - ring_atoms
        rn = len(set(rings))
        if rn == 0:
            return f"C{heavy_atoms}00"
        elif rn == 1:
            rp = ring_p[0]
            if set(large_flag_arr) == {True}:
                return f"RA{__arom_number[rp]}{out_ring_atoms}"
            else:
                return f"R{ring_atoms}{__arom_number[rp]}{out_ring_atoms}"
        else:
            return f"RF{rn}{out_ring_atoms}"

    def elementary_frags(self, m):
        elem_frag = {}
        for i in range(len(self.frags)):
            elem_frag[i] = {
                "index": i,
                "components": self.frags[i],
                "offsets": [],
                "connects": [],
                "label": self.elem_frag_label(m, self.frags[i]),
            }
            for a in self.frags[i]:
                n = sum(m.Atoms[a].break_bond)
                if n > 0:
                    elem_frag[i]["offsets"].append([a, n, "H"])
                    ind = [ii for ii, x in enumerate(m.Atoms[a].break_bond) if x]
                    for j in ind:
                        b = m.Atoms[a].connect[j]
                        for jj in range(len(self.frags)):
                            if b in self.frags[jj]:
                                elem_frag[i]["connects"].append([jj, a, b])
        for index, frag in elem_frag.items():
            if frag["label"][0] == "C":
                r_frag = [
                    rc[0]
                    for rc in frag["connects"]
                    if elem_frag[rc[0]]["label"][0] == "R"
                ]
                if len(r_frag) == 0:
                    frag["ef_type"] = "chain"
                elif len(r_frag) == 1:
                    frag["ef_type"] = "side"
                else:
                    frag["ef_type"] = "link"
            else:
                frag["ef_type"] = "ring"
        m.elem_frag = elem_frag

    def sssr_frags(self, m):
        __arom_number = {
            "nonar": 0,
            "ar1": 1,
            "ar2": 2,
            "ar3": 3,
            "ar4": 4,
            "ar5": 5,
            "ar6": 6,
            "ar7": 7,
            "ar8": 8,
        }
        sssr_frag = {}
        ii = -1
        for __, ring in m.ring_dict.items():
            ii += 1
            sssr_frag[ii] = {
                "index": ii,
                "components": ring[:-1],
                "offsets": [],
                "special_bond": [],
            }
            sssr_frag[ii]["label"] = (
                f"R{len(ring) - 1}{__arom_number[ring[-1]]}0"
                if len(ring) - 1 < 10
                else f"RA{__arom_number[ring[-1]]}0"
            )
            for a in ring[:-1]:
                double_in_other_ring = []
                n = 0
                atom_a = m.Atoms[a]
                for j in range(len(atom_a.connect)):
                    atom_b = m.Atoms[atom_a.connect[j]]
                    if atom_b.No not in sssr_frag[ii]["components"]:
                        if not atom_a.break_bond[j]:
                            if len(atom_b.inring) == 0:
                                sssr_frag[ii]["components"].append(atom_b.No)
                            else:
                                n += 1
                                if atom_a.bond_type[j] == "2":
                                    double_in_other_ring.append(atom_b.No)
                        else:
                            n += 1
                if n > 0:
                    sssr_frag[ii]["offsets"].append([a, n, "H"])
                if len(double_in_other_ring) > 0:
                    sssr_frag[ii]["special_bond"].append(a)
        m.sssr_frag = sssr_frag

    def get_ring_renumber(self, dicts, cc, ii, tmp):
        aa = list(set(dicts[cc]).intersection(set(dicts.keys())).difference(set(tmp)))
        if len(aa) > 0:
            ii += 1
            tmp.append(aa[0])
            self.renumber_dict_tmp[aa[0]] = ii
            self.get_ring_renumber(dicts, aa[0], ii, tmp)

    def get_chain_renumber(self, m, frag):
        atoms = list(
            set(frag["components"]).difference(
                set(
                    [
                        aa
                        for aa in frag["components"]
                        if m.Atoms[aa].elem in ["H", "F", "Cl", "Br", "I"]
                        and len(m.Atoms[aa].connect) == 1
                    ]
                )
            )
        )
        connect_dict = {}
        for aa in atoms:
            connect_dict[aa] = list(set(m.Atoms[aa].connect).intersection(set(atoms)))
        terminal_atom_arr = [aa for aa, bb in connect_dict.items() if len(bb) == 1]
        has_ring_dict = {}
        for aa in frag["components"]:
            has_ring_dict[aa] = [
                rc
                for rc in frag["connects"]
                if rc[1] == aa and m.elem_frag[rc[0]]["label"][0] == "R"
            ]
        for aa, bb in has_ring_dict.items():
            has_ring_dict[aa] = len(bb)
        if len(terminal_atom_arr) == 0:
            terminal_atom_arr = [atoms[0], atoms[1]]
        elif len(terminal_atom_arr) == 1:
            for aa in atoms:
                if aa != terminal_atom_arr[0]:
                    terminal_atom_arr.append(aa)
                    break
        all_chain = search_chain(
            terminal_atom_arr, connect_dict, additional_conditions=has_ring_dict
        )
        renumber_dict = {}
        ii = -1
        for aa in all_chain[0]:
            ii += 1
            renumber_dict[aa] = ii
        res_atoms = set(frag["components"]).difference(set(renumber_dict.keys()))
        for aa in res_atoms:
            tmp = list(set(m.Atoms[aa].connect).intersection(set(renumber_dict.keys())))
            if len(tmp) > 0:
                renumber_dict[aa] = renumber_dict[tmp[0]]
        res_atoms = set(frag["components"]).difference(set(renumber_dict.keys()))
        for aa in res_atoms:
            renumber_dict[aa] = "X"
        return renumber_dict

    def get_renumber(self, m, frag):
        renumber_dict = {}
        heavy_atom = list(
            set(frag["components"]).difference(
                set(
                    [
                        aa
                        for aa in frag["components"]
                        if m.Atoms[aa].elem in ["H", "F", "Cl", "Br", "I"]
                        and len(m.Atoms[aa].connect) == 1
                    ]
                )
            )
        )
        if frag["label"][0] == "C":
            ER_atoms = [
                cc
                for cc in frag["components"]
                if m.Atoms[cc].local in ["LMC", "LTC", "LT"]
            ]
            if len(ER_atoms) == 1 or len(heavy_atom) == 1:
                for cc in frag["components"]:
                    renumber_dict[cc] = 0
                return renumber_dict
            if len(heavy_atom) == 2:
                renumber_dict[heavy_atom[0]] = 0
                renumber_dict[heavy_atom[1]] = 1
                return renumber_dict
            if len(heavy_atom) >= 10:
                for cc in frag["components"]:
                    renumber_dict[cc] = "X"
                return renumber_dict
            return self.get_chain_renumber(m, frag)
        elif frag["label"][0] == "R" and frag["label"][1] == "F":
            for cc in frag["components"]:
                renumber_dict[cc] = "X"
            return renumber_dict
        else:
            out_ring = {}
            in_ring = {}
            cb = 0
            for ca in heavy_atom:
                if len(m.Atoms[ca].inring) > 0:
                    in_ring[ca] = m.Atoms[ca].connect
                    cb = ca
            for ca in heavy_atom:
                if len(m.Atoms[ca].inring) == 0:
                    try:
                        out_ring[ca] = list(
                            set(m.Atoms[ca].connect).intersection(set(in_ring.keys()))
                        )[0]
                    except:
                        out_ring[ca] = list(
                            set(m.Atoms[ca].connect).intersection(set(heavy_atom))
                        )[0]
            self.renumber_dict_tmp = {}
            self.renumber_dict_tmp[cb] = 0
            self.get_ring_renumber(in_ring, cb, 0, [cb])
            renumber_dict = self.renumber_dict_tmp
            for out_atom, in_atom in out_ring.items():
                try:
                    renumber_dict[out_atom] = renumber_dict[in_atom]
                except:
                    renumber_dict[out_atom] = "X"
            return renumber_dict

    def get_offset_nonef(self, frags):
        offsets = []
        indexs = []
        for frag in frags:
            if isinstance(frag["index"], list):
                indexs.extend(frag["index"])
            else:
                indexs.append(frag["index"])
        total_offset = []
        for frag in frags:
            total_offset.extend(frag["offsets"])
        del_offset_node = {}
        for frag in frags:
            for rr in frag["connects"]:
                if rr[0] in indexs:
                    if rr[1] not in del_offset_node:
                        del_offset_node[rr[1]] = 0
                    del_offset_node[rr[1]] += 1
        for of in total_offset:
            if of[0] in del_offset_node:
                nn = of[1] - del_offset_node[of[0]]
                if nn < 1:
                    pass
                else:
                    offsets.append([of[0], nn, of[2]])
            else:
                offsets.append(of)
        return offsets
        for frag in frags:
            offsets += [rr for rr in frag["offsets"] if rr[0] not in indexs]
        return offsets

    def create_null_frag(self, indexs):
        frag_dict = {
            "components": [],
            "offsets": [],
            "label": "C000",
            "connects": [],
            "index": -1,
        }
        for ii in indexs:
            frag_dict["connects"].append([ii, -1, -1])
        return frag_dict

    def get_nonef_info(self, m, frags, frag_type, extend_frags=None):
        if extend_frags is None:
            extend_frags = []
        info_dict = {"components": [], "connects": []}
        if frag_type != "tf":
            info_dict["label"] = self.get_nonef_frag_label(m, frags)
        else:
            info_dict["label"] = frags[0]["tmp_label"]
        for frag in frags + extend_frags:
            info_dict["components"].extend(frag["components"])
            info_dict["connects"].extend(frag["connects"])
        info_dict["offsets"] = self.get_offset_nonef(frags + extend_frags)
        return info_dict

    def get_nonef_frag_label(self, m, frags):
        frag_type_arr = [frag["label"][0] for frag in frags]
        if len(set(frag_type_arr)) == 1 and "C" in frag_type_arr:
            nn = 0
            for frag in frags:
                nn += int(frag["label"][1])
            if nn >= 10:
                nn = "A"
            return f"C{nn}00"
        if frags[0]["label"][0] == "C" and frags[0]["label"][1] == "0":
            label = f"C00{len(frags[1:])}_"
            for frag in frags[1:]:
                label += frag["label"][:3] + "0"
            return label
        renumber_dict = self.get_renumber(m, frags[0])
        label = ""
        nn = int(frags[0]["label"][3]) + len(frags[1:])
        label = (
            frags[0]["label"][:3] + str(nn) if nn < 10 else frags[0]["label"][:3] + "A"
        )
        label += "_"
        for frag in frags[1:]:
            connect_atom = [
                rr[1] for rr in frags[0]["connects"] if rr[0] == frag["index"]
            ]
            label += (
                frag["label"][:3] + str(renumber_dict[connect_atom[0]])
                if "tmp_label" not in frag.keys()
                else frag["tmp_label"][:3] + str(renumber_dict[connect_atom[0]])
            )
        return label

    def csf_frags(self, m, tmp_frag_index=None):
        if tmp_frag_index is None:
            tmp_frag_index = []
        csf_frag = {}
        if len(m.elem_frag) == 2 and any(
            connection[0] == 1 for connection in m.elem_frag[0]["connects"]
        ):
            tmp_frag_index.append(set([0, 1]))
            csf_frag[0] = {
                "components": [i for i in range(len(m.Atoms))],
                "offsets": [],
                "label": self.get_nonef_frag_label(m, [m.elem_frag[0], m.elem_frag[1]]),
            }
        else:
            ii = -1
            for index, frag in m.elem_frag.items():
                if frag["label"][0] == "C":
                    c_frag = [
                        rr[0]
                        for rr in frag["connects"]
                        if m.elem_frag[rr[0]]["label"][0] == "C"
                    ]
                    r_frag = [
                        rr[0]
                        for rr in frag["connects"]
                        if m.elem_frag[rr[0]]["label"][0] == "R"
                    ]
                    for rc in c_frag:
                        if set([index, rc]) not in tmp_frag_index:
                            tmp_frag_index.append(set([index, rc]))
                            ii += 1
                            csf_frag[ii] = self.get_nonef_info(
                                m, [frag, m.elem_frag[rc]], "csf"
                            )
                    if len(r_frag) > 0:
                        frag_indexs = [index]
                        frag_indexs.extend([rc for rc in r_frag])
                        if set(frag_indexs) not in tmp_frag_index:
                            frags = [frag]
                            frags.extend([m.elem_frag[rc] for rc in r_frag])
                            tmp_frag_index.append(set(frag_indexs))
                            ii += 1
                            csf_frag[ii] = self.get_nonef_info(m, frags, "csf")
                else:
                    for rr in frag["connects"]:
                        if m.elem_frag[rr[0]]["label"][0] == "R":
                            if set([index, rr[0]]) not in tmp_frag_index:
                                tmp_frag_index.append(set([index, rr[0]]))
                                ii += 1
                                null_frag = self.create_null_frag([rr[0]])
                                csf_frag[ii] = self.get_nonef_info(
                                    m, [null_frag, frag, m.elem_frag[rr[0]]], "csf"
                                )
        m.csf_frag = csf_frag
        return tmp_frag_index

    def rsf_frags(self, m, tmp_frag_index=None, long_ring=1, *, build_products=True):
        if tmp_frag_index is None:
            tmp_frag_index = []
        ii = -1
        rsf_frag = {}
        for index, frag in m.elem_frag.items():
            if frag["label"][0] == "R":
                frag_indexs = [index]
                frag_indexs.extend([rr[0] for rr in frag["connects"]])
                extend_indexs = []
                for rr in frag["connects"]:
                    if (
                        m.elem_frag[rr[0]]["label"][0] == "C"
                        and int(m.elem_frag[rr[0]]["label"][1]) <= long_ring
                    ):
                        tmp_frag = m.elem_frag[rr[0]]
                        tmp_connect = [
                            rs[0]
                            for rs in tmp_frag["connects"]
                            if m.elem_frag[rs[0]]["label"][0] == "R" and rs[0] != index
                        ]
                        if len(tmp_connect) > 0:
                            m.elem_frag[rr[0]]["tmp_label"] = (
                                m.elem_frag[rr[0]]["label"][:2] + "R0"
                            )
                            extend_indexs.extend(tmp_connect)
                if build_products and len(frag_indexs) > 2:
                    if set(frag_indexs + extend_indexs) not in tmp_frag_index:
                        tmp_frag_index.append(set(frag_indexs))
                        frags = [m.elem_frag[rc] for rc in frag_indexs]
                        extend_frags = [m.elem_frag[rc] for rc in set(extend_indexs)]
                        ii += 1
                        rsf_frag[ii] = self.get_nonef_info(
                            m, frags, "rsf", extend_frags=extend_frags
                        )
        m.rsf_frag = rsf_frag
        return tmp_frag_index

    def chain_frags(self, cc, m):
        c_frag = set(
            [
                rr[0]
                for rr in m.elem_frag[cc]["connects"]
                if m.elem_frag[rr[0]]["label"][0] == "C"
            ]
        ).difference(set(self.chain_frag_tmp))
        for index in c_frag:
            self.chain_frag_tmp.append(index)
            self.chain_frags(index, m)

    def scaffold_tf_frags(self, m, *, include_scaffold=True, include_tf=True):
        if include_tf and not include_scaffold:
            raise ValueError("TF construction requires scaffold intermediates")
        scaffold_frag = {}
        chain_frag = {}
        tf_frag = {}
        cef_used = []
        cii = -1
        sii = -1
        tii = -1
        for index, frag in m.elem_frag.items():
            if frag["label"][0] == "C" and index not in cef_used:
                self.chain_frag_tmp = [index]
                self.chain_frags(index, m)
                cef_used.extend(self.chain_frag_tmp)
                cii += 1
                frags = [m.elem_frag[index] for index in self.chain_frag_tmp]
                chain_frag[cii] = self.get_nonef_info(m, frags, "chainf")
                chain_frag[cii]["index"] = self.chain_frag_tmp
                chain_frag[cii]["connects"] = []
                for frag in frags:
                    chain_frag[cii]["connects"].extend(
                        [
                            rr
                            for rr in frag["connects"]
                            if rr[0] not in self.chain_frag_tmp
                        ]
                    )
                if not include_scaffold:
                    continue
                frags = [chain_frag[cii]]
                tmp_indexs = []
                for rc in self.chain_frag_tmp:
                    tmp_indexs.extend(
                        [
                            rs[0]
                            for rs in m.elem_frag[rc]["connects"]
                            if m.elem_frag[rs[0]]["label"][0] == "R"
                        ]
                    )
                    frags.extend(
                        [
                            m.elem_frag[rs[0]]
                            for rs in m.elem_frag[rc]["connects"]
                            if m.elem_frag[rs[0]]["label"][0] == "R"
                        ]
                    )
                if len(frags) > 2:
                    sii += 1
                    scaffold_frag[sii] = self.get_nonef_info(m, frags, "scaffoldf")
                    scaffold_frag[sii]["index"] = self.chain_frag_tmp + tmp_indexs
                    scaffold_frag[sii]["tmp_label"] = scaffold_frag[sii]["label"] + "_"
                    scaffold_frag[sii]["connects"] = []
                    for frag in frags[1:]:
                        scaffold_frag[sii]["connects"].extend(
                            [
                                rr
                                for rr in frag["connects"]
                                if rr[0] not in self.chain_frag_tmp
                            ]
                        )
                    extend_frags = []
                    for frag in frags[1:]:
                        tmp_frags = [
                            m.elem_frag[rc[0]]
                            for rc in frag["connects"]
                            if rc[0] not in self.chain_frag_tmp
                        ]
                        scaffold_frag[sii]["tmp_label"] += str(len(tmp_frags))
                        extend_frags.extend(tmp_frags)
                    if include_tf:
                        tii += 1
                        tf_frag[tii] = self.get_nonef_info(
                            m, [scaffold_frag[sii]], "tf", extend_frags=extend_frags
                        )
            elif include_scaffold and frag["label"][0] == "R" and index not in cef_used:
                cef_used.append(index)
                for rr in frag["connects"]:
                    if m.elem_frag[rr[0]]["label"][0] == "R":
                        if rr[0] not in cef_used:
                            cef_used.append(rr[0])
                            sii += 1
                            null_frag = self.create_null_frag([rr[0]])
                            scaffold_frag[sii] = self.get_nonef_info(
                                m, [null_frag, frag, m.elem_frag[rr[0]]], "scaffoldf"
                            )
                            scaffold_frag[sii]["index"] = [index, rr[0]]
                            scaffold_frag[sii]["tmp_label"] = (
                                scaffold_frag[sii]["label"] + "_"
                            )
                            extend_frags = []
                            tmp_frags = [
                                m.elem_frag[rc[0]]
                                for rc in frag["connects"]
                                if rc[0] != rr[0]
                            ]
                            scaffold_frag[sii]["tmp_label"] += str(len(tmp_frags))
                            extend_frags.extend(tmp_frags)
                            tmp_frags = [
                                m.elem_frag[rc[0]]
                                for rc in m.elem_frag[rr[0]]["connects"]
                                if rc[0] != index
                            ]
                            scaffold_frag[sii]["tmp_label"] += str(len(tmp_frags))
                            extend_frags.extend(tmp_frags)
                            if include_tf:
                                tii += 1
                                tf_frag[tii] = self.get_nonef_info(
                                    m,
                                    [scaffold_frag[sii]],
                                    "tf",
                                    extend_frags=extend_frags,
                                )
            else:
                cef_used.append(index)
        m.chain_frag = chain_frag
        m.scaffold_frag = scaffold_frag
        m.tf_frag = tf_frag

    def secondary_frags(self, m):
        seco_frag = {}
        if len(m.elem_frag) == 2 and any(
            connection[0] == 1 for connection in m.elem_frag[0]["connects"]
        ):
            seco_frag[0] = {
                "components": [i for i in range(len(m.Atoms))],
                "offsets": [],
            }
        else:
            ii = -1
            for index, frag in m.elem_frag.items():
                n = len(frag["connects"])
                if n > 1:
                    ii += 1
                    seco_frag[ii] = {}
                    seco_frag[ii]["components"] = deepcopy(frag["components"])
                    seco_frag[ii]["offsets"] = deepcopy(frag["offsets"])
                    combine_a = []
                    for cc in frag["connects"]:
                        fn = cc[0]
                        combine_a.append(cc[1])
                        combine_a.append(cc[2])
                        seco_frag[ii]["components"] += m.elem_frag[fn]["components"]
                        seco_frag[ii]["offsets"] += m.elem_frag[fn]["offsets"]
                    for a in set(combine_a):
                        for iii in range(len(seco_frag[ii]["offsets"])):
                            if seco_frag[ii]["offsets"][iii][0] == a:
                                break
                        del seco_frag[ii]["offsets"][iii]
        m.seco_frag = seco_frag

    def sketch_frags(self, m):
        sketch_frag = {}
        nn = max([index for index in m.chain_frag]) if len(m.chain_frag) > 0 else 0
        for index, frag in m.elem_frag.items():
            if frag["ef_type"] == "ring":
                nn += 1
                sketch_frag[nn] = deepcopy(frag)
                sketch_frag[nn]["index"] = [index]
                sketch_frag[nn]["type"] = "Ring"
                sketch_frag[nn]["frag_connects"] = []
        for index, frag in m.chain_frag.items():
            sketch_frag[index] = deepcopy(frag)
            if len(frag["connects"]) > 1:
                sketch_frag[index]["type"] = "Link"
            else:
                sketch_frag[index]["type"] = "Side"
            sketch_frag[index]["frag_connects"] = []
            for rr in frag["connects"]:
                r_index = [
                    index
                    for index, frag in sketch_frag.items()
                    if rr[0] in frag["index"]
                ][0]
                ra = deepcopy(rr)
                ra[0] = r_index
                sketch_frag[index]["frag_connects"].append(ra)
        for index, frag in sketch_frag.items():
            if frag["type"] == "Ring":
                for rr in frag["connects"]:
                    r_index = [
                        index
                        for index, frag in sketch_frag.items()
                        if rr[0] in frag["index"]
                    ][0]
                    ra = deepcopy(rr)
                    ra[0] = r_index
                    sketch_frag[index]["frag_connects"].append(ra)
        m.sketch_frag = sketch_frag
