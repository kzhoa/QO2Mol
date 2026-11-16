import math

import numpy as np

import logging
logger = logging.getLogger(__name__)
from ....graph import Molecule
from ....shared.geometry import Angle, Dihedral


class ZMatrix:
    def __init__(self):
        self.atoms = []  # extra attr: pos_array
        self.partners = []
        self.bonds = [None]  # extra attr: type
        self.angles = [None, None]  # extra attr: type21, type23
        self.dihedrals = [None, None, None]

    def from_molecule(self, mol: Molecule):
        self.atoms = mol.Atoms[:]
        self.partners = [[] for _ in self.atoms]

        # remove ring-closing bonds and ensure a1 < a2
        for a2 in range(1, len(self.atoms)):
            try:
                a1 = max(filter(lambda x: x < a2, mol.connectivity[a2]))
            except ValueError:
                logger.error("Molecule cannot be converted to ZMatrix preserving atom order")
                raise

            bond = next(bond for bond in mol.Bonds if {a1, a2} == {bond.a1, bond.a2})
            bond.calc_value([mol.Atoms[i].coor for i in (a1, a2)])
            bond.type = mol.Atoms[a2].bond_type[mol.connectivity[a2].index(a1)]
            self.bonds.append(bond)

        # new connectivity without ring-closing bonds
        for bond in self.bonds[1:]:
            self.partners[bond.a1].append(bond.a2)
            self.partners[bond.a2].append(bond.a1)

        for a3 in range(2, len(self.atoms)):
            bond23 = self.bonds[a3]
            a2 = bond23.a1
            a1 = min(self.partners[a2])

            angle = Angle("harmonic", a1, a2, a3)
            angle.calc_value([mol.Atoms[i].coor for i in (a1, a2, a3)])
            angle.type21 = mol.Atoms[a2].bond_type[mol.connectivity[a2].index(a1)]
            angle.type23 = mol.Atoms[a2].bond_type[mol.connectivity[a2].index(a3)]
            self.angles.append(angle)

        for a4 in range(3, len(self.atoms)):
            angle234 = self.angles[a4]
            a2, a3 = angle234.a1, angle234.a2
            partners = [p for p in self.partners[a2] if p not in (a3, a4)]
            if partners == [] or min(partners) > a4:
                partners = [p for p in self.partners[a3] if p not in (a2, a3)]
            a1 = min(partners)
            dihedral = Dihedral("amber", a1, a2, a3, a4)
            dihedral.calc_value([mol.Atoms[i].coor for i in (a1, a2, a3, a4)])
            self.dihedrals.append(dihedral)

        for atom in self.atoms:
            atom.pos_array = np.array(atom.coor, dtype=np.float64)

    def to_str(self):
        string = ""
        for atom, bond, angle, dihedral in zip(self.atoms, self.bonds, self.angles, self.dihedrals):
            string += f"{atom.elem:2s}"
            if bond is not None:
                string += f" {(bond.a1 + 1):4d} {bond.value:10.4f}"
            if angle is not None:
                string += f" {(angle.a1 + 1):4d} {angle.value_a:10.4f}"
            if dihedral is not None:
                string += f" {(dihedral.a1 + 1):4d} {dihedral.value_a:10.4f}"
            string += "\n"

        return string

    def to_coordinates(self):
        if len(self.atoms) >= 1:
            self.atoms[0].pos_array[:] = [0, 0, 0]

        if len(self.atoms) >= 2:
            self.atoms[1].pos_array[:] = [self.bonds[1].value, 0, 0]

        # third atom at distance r from ir forms angle a 3-ir-ia in plane xy
        if len(self.atoms) >= 3:
            r = self.bonds[2].value
            ir = self.bonds[2].a1
            ang = self.angles[2].value
            ia = self.angles[2].a1

            # for this construction, the new atom is at point (x, y), atom
            # ir is at point (xr, yr) and atom ia is at point (xa, ya).
            # Theta is the angle between the vector joining ir to ia and
            # the x-axis, a' (= theta - a) is is the angle between r and
            # the x-axis. x = xa + r cos a', y = ya + r sin a'.  From the
            # dot product of a unitary vector along x with the vector from
            # ir to ia, theta can be calculated: cos theta = (xa - xr) /
            # sqrt((xa - xr)^2 + (ya - yr)^2).  If atom ia is in third or
            # forth quadrant relative to atom ir, ya - yr < 0, then theta
            # = 2 pi - theta. */
            dx, dy, dz = self.atoms[ia].pos_array - self.atoms[ir].pos_array
            theta = math.acos(np.clip(dx / math.sqrt(dx * dx + dy * dy), -1, 1))
            if dy < 0.0:
                theta = 2 * math.pi - theta
            ang = theta - ang
            self.atoms[2].pos_array[:] = [
                self.atoms[ir].pos_array[0] + r * math.cos(ang),
                self.atoms[ir].pos_array[1] + r * math.sin(ang),
                0.0,
            ]

        # nth atom at distance r from atom ir forms angle a at 3-ir-ia
        # and dihedral angle between planes 3-ir-ia and ir-ia-id
        if len(self.atoms) >= 4:
            for i in range(3, len(self.atoms)):
                r = self.bonds[i].value
                ir = self.bonds[i].a1
                ang = self.angles[i].value
                ia = self.angles[i].a1
                dih = self.dihedrals[i].value
                id = self.dihedrals[i].a1

                # for this construction the new atom is at point A, atom ir is
                # at B, atom ia at C and atom id at D.  Point a is the
                # projection of A onto the plane BCD.  Point b is the
                # projection of A along the direction BC (the line defining
                # the dihedral angle between planes ABC and BCD). n = CD x BC
                # / |CD x BC| is the unit vector normal to the plane BCD. m =
                # BC x n / |BC x n| is the unit vector on the plane BCD normal
                # to the direction BC.
                #
                #                               .'A
                #                 ------------.' /.-----------------
                #                /           b /  .               /
                #               /           ./    .              /
                #              /           B......a      ^      /
                #             /           /              |n    /
                #            /           /                    /
                #           /           C                    /
                #          /             \                  /
                #         /               \                /
                #        /plane BCD        D              /
                #       ----------------------------------
                #
                #                    A              C------B...b
                #                   /.             /        .  .
                #                  / .            /    |m    . .
                #                 /  .           /     V      ..
                #         C------B...b          D              a
                #

                BA = r
                vB = self.atoms[ir].pos_array
                vC = self.atoms[ia].pos_array
                vD = self.atoms[id].pos_array

                vBC = vC - vB
                vCD = vD - vC

                BC = math.sqrt(vBC.dot(vBC))
                bB = BA * math.cos(ang)
                bA = BA * math.sin(ang)
                aA = bA * math.sin(dih)
                ba = bA * math.cos(dih)

                vb = vC - vBC * ((BC - bB) / BC)
                vn = np.cross(vCD, vBC)
                vn = vn / math.sqrt(vn.dot(vn))
                vm = np.cross(vBC, vn)
                vm = vm / math.sqrt(vm.dot(vm))
                va = vb + vm * ba
                vA = va + vn * aA

                self.atoms[i].pos_array = vA

        return np.array([atom.pos_array for atom in self.atoms])
