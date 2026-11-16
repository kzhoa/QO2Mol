from dataclasses import dataclass
from typing import List, Union

class Freezeable:
    def _freeze(self):
        self.__is_frozen = True

    def __setattr__(self, name: str, value) -> None:
        if not hasattr(self, "_Freezeable__is_frozen"):
            object.__setattr__(self, "_Freezeable__is_frozen", False)
        if name != "_Freezeable__is_frozen" and self.__is_frozen and not hasattr(self, name):
            raise KeyError("Cannot create new attribute in frozen object")
        object.__setattr__(self, name, value)


class ElementError(ValueError):
    def __init__(self, element):
        super().__init__(f"Element {element} is not a legal element.")

@dataclass
class Element(Freezeable):
    number: int  # 编号
    symbol: str  # 符号
    mass: float  # 原子量
    valent_radius: float # 共价半径(单键)
    double_bond_radius: float = None # 元素形成双键的半径
    triple_bond_radius: float = None # 元素形成三键的半径
    conju_bond_radius: float = None # 元素形成共轭键的半径
    vdw_radius: float = None  # 范德华半径
    valents: List[int] = None  # 化合价

    def __post_init__(self):
        if self.valents is None:
            self.valents = []
        self._freeze()

    @staticmethod
    def get(key: Union[float, int, str], check_mass=False):
        if key not in _element_table:  # From number or symbol
            if check_mass and (type(key) in [int, float]):
                return Element.get_by_mass(key)
            raise ElementError(key)
        return _element_table[key]

    @staticmethod
    def get_by_mass(mass: Union[float, int]):
        for element in ELEMENT_TABLE:
            if abs(mass - element.mass) <= 0.1:
                return element
        else:
            raise ElementError(f"with mass {mass}")

    def __str__(self) -> str:
        return f"<{self.symbol} {self.number}>"


# radius_vdw from https://pubs.acs.org/doi/10.1021/jp8111556
# valent_radius from https://periodictable.com/Properties/A/CovalentRadius.an.html

ELEMENT_TABLE: List[Element] = [
    Element(number=1, symbol="H", mass=1.008, valent_radius=0.3, vdw_radius=1.10,valents=[1]),
    Element(number=2, symbol="He", mass=4.003, valent_radius=1.16,valents=[0]),
    Element(number=3, symbol="Li", mass=6.941, valent_radius=1.23, valents=[1]),
    Element(number=4, symbol="Be", mass=9.012, valent_radius=0.89, valents=[2]),
    Element(number=5, symbol="B", mass=10.811, valent_radius=0.88, vdw_radius=1.92,valents=[3]),
    Element(number=6, symbol="C", mass=12.011, valent_radius=0.77, vdw_radius=1.70,
            double_bond_radius=0.6776, conju_bond_radius=0.7007,triple_bond_radius=0.6006,valents=[4]),
    Element(number=7, symbol="N", mass=14.007, valent_radius=0.7, vdw_radius=1.55,
            double_bond_radius=0.616, conju_bond_radius=0.637, triple_bond_radius=0.546,valents=[3, 5]),
    Element(number=8, symbol="O", mass=15.999, valent_radius=0.66, vdw_radius=1.52,
            double_bond_radius=0.5808,conju_bond_radius=0.6006,valents=[2]),
    Element(number=9, symbol="F", mass=18.988, valent_radius=0.58, vdw_radius=1.47,valents=[1]),
    Element(number=10, symbol="Ne", mass=20.17, valent_radius=0.55),
    Element(number=11, symbol="Na", mass=22.99, valent_radius=1.4, valents=[1]),
    Element(number=12, symbol="Mg", mass=24.305, valent_radius=1.36, valents=[2]),
    Element(number=13, symbol="Al", mass=26.982, valent_radius=1.25, valents=[3, 5]),
    Element(number=14, symbol="Si", mass=28.085, valent_radius=1.17, vdw_radius=2.07,
            double_bond_radius=1.0296, triple_bond_radius=0.9126,valents=[4]),
    Element(number=15, symbol="P", mass=30.974, valent_radius=1.05, vdw_radius=1.80,
            double_bond_radius=0.924,conju_bond_radius=0.9555,valents=[3, 5, 4]),
    Element(number=16, symbol="S", mass=32.06, valent_radius=1.01, vdw_radius=1.80,
            double_bond_radius=0.9624,valents=[2, 4, 6]), ##开始的值double_bond_radius= 0.92
    Element(number=17, symbol="Cl", mass=35.453, valent_radius=0.99, vdw_radius=1.75, valents=[1]),
    Element(number=18, symbol="Ar", mass=39.94, valent_radius=1.55),
    Element(number=19, symbol="K", mass=39.089, valent_radius=2.03, valents=[1]),
    Element(number=20, symbol="Ca", mass=40.08, valent_radius=1.74, valents=[2]),
    Element(number=21, symbol="Sc", mass=44.956, valent_radius=1.44),
    Element(number=22, symbol="Ti", mass=47.9, valent_radius=1.32),
    Element(number=23, symbol="V", mass=50.941, valent_radius=1.2),
    Element(number=24, symbol="Cr", mass=51.996, valent_radius=1.13),
    Element(number=25, symbol="Mn", mass=54.938, valent_radius=1.17),
    Element(number=26, symbol="Fe", mass=55.84, valent_radius=1.16),
    Element(number=27, symbol="Co", mass=58.933, valent_radius=1.16),
    Element(number=28, symbol="Ni", mass=58.69, valent_radius=1.15),
    Element(number=29, symbol="Cu", mass=63.54, valent_radius=1.17),
    Element(number=30, symbol="Zn", mass=65.38, valent_radius=1.25),
    Element(number=31, symbol="Ga", mass=69.72, valent_radius=1.25),
    Element(number=32, symbol="Ge", mass=72.59, valent_radius=1.22),
    Element(number=33, symbol="As", mass=74.922, valent_radius=1.21, vdw_radius=2.05),
    Element(number=34, symbol="Se", mass=78.9, valent_radius=1.17, vdw_radius=2.29),
    Element(number=35, symbol="Br", mass=79.904, valent_radius=1.20, vdw_radius=1.83, valents=[1]),
    Element(number=36, symbol="Kr", mass=83.8, valent_radius=1.89),
    Element(number=37, symbol="Rb", mass=85.467, valent_radius=2.25),
    Element(number=38, symbol="Sr", mass=87.62, valent_radius=1.92),
    Element(number=39, symbol="Y", mass=88.906, valent_radius=1.62),
    Element(number=40, symbol="Zr", mass=91.22, valent_radius=1.45),
    Element(number=41, symbol="Nb", mass=92.906, valent_radius=1.34),
    Element(number=42, symbol="Mo", mass=95.94, valent_radius=1.29),
    Element(number=43, symbol="Tc", mass=99.0, valent_radius=1.23),
    Element(number=44, symbol="Ru", mass=101.07, valent_radius=1.24),
    Element(number=45, symbol="Rh", mass=102.906, valent_radius=1.25),
    Element(number=46, symbol="Pd", mass=106.42, valent_radius=1.28),
    Element(number=47, symbol="Ag", mass=107.868, valent_radius=1.34),
    Element(number=48, symbol="Cd", mass=112.41, valent_radius=1.41),
    Element(number=49, symbol="In", mass=114.82, valent_radius=1.5),
    Element(number=50, symbol="Sn", mass=118.6, valent_radius=1.4),
    Element(number=51, symbol="Sb", mass=121.7, valent_radius=1.41),
    Element(number=52, symbol="Te", mass=127.6, valent_radius=1.37, vdw_radius=2.5),
    Element(number=53, symbol="I", mass=126.905, valent_radius=1.39, vdw_radius=1.98, valents=[1]),
    Element(number=54, symbol="Xe", mass=131.3, valent_radius=2.09),
    Element(number=55, symbol="Cs", mass=132.905, valent_radius=2.35),
    Element(number=56, symbol="Ba", mass=137.33, valent_radius=1.98),
    Element(number=57, symbol="La", mass=138.905, valent_radius=1.69),
    Element(number=58, symbol="Ce", mass=140.12, valent_radius=1.65),
    Element(number=59, symbol="Pr", mass=140.91, valent_radius=1.65),
    Element(number=60, symbol="Nd", mass=144.2, valent_radius=1.64),
    Element(number=61, symbol="Pm", mass=147.0, valent_radius=1.64),
    Element(number=62, symbol="Sm", mass=150.4, valent_radius=1.66),
    Element(number=63, symbol="Eu", mass=151.96, valent_radius=1.85),
    Element(number=64, symbol="Gd", mass=157.25, valent_radius=1.61),
    Element(number=65, symbol="Tb", mass=158.93, valent_radius=1.59),
    Element(number=66, symbol="Dy", mass=162.5, valent_radius=1.59),
    Element(number=67, symbol="Ho", mass=164.93, valent_radius=1.58),
    Element(number=68, symbol="Er", mass=167.2, valent_radius=1.57),
    Element(number=69, symbol="Tm", mass=168.934, valent_radius=1.56),
    Element(number=70, symbol="Yb", mass=173.0, valent_radius=1.7),
    Element(number=71, symbol="Lu", mass=174.96, valent_radius=1.56),
    Element(number=72, symbol="Hf", mass=178.4, valent_radius=1.44),
    Element(number=73, symbol="Ta", mass=180.947, valent_radius=1.34),
    Element(number=74, symbol="W", mass=183.8, valent_radius=1.3),
    Element(number=75, symbol="Re", mass=186.207, valent_radius=1.28),
    Element(number=76, symbol="Os", mass=190.2, valent_radius=1.26),
    Element(number=77, symbol="Ir", mass=192.2, valent_radius=1.26),
    Element(number=78, symbol="Pt", mass=195.08, valent_radius=1.29),
    Element(number=79, symbol="Au", mass=196.967, valent_radius=1.34),
    Element(number=80, symbol="Hg", mass=200.5, valent_radius=1.44),
    Element(number=81, symbol="Tl", mass=204.3, valent_radius=1.55),
    Element(number=82, symbol="Pb", mass=207.2, valent_radius=1.54),
    Element(number=83, symbol="Bi", mass=208.98, valent_radius=1.52),
    Element(number=84, symbol="Po", mass=209.0, valent_radius=1.53),
    Element(number=85, symbol="At", mass=210.0, valent_radius=1.52),
    Element(number=86, symbol="Rn", mass=222.0, valent_radius=1.53),
    Element(number=87, symbol="Fr", mass=223.0, valent_radius=2.45),
    Element(number=88, symbol="Ra", mass=226.03, valent_radius=2.02),
    Element(number=89, symbol="Ac", mass=227.0, valent_radius=1.7),
    Element(number=90, symbol="Th", mass=232.03, valent_radius=1.63),
    Element(number=91, symbol="Pa", mass=231.03, valent_radius=1.46),
    Element(number=92, symbol="U", mass=238.02, valent_radius=1.4),
    Element(number=0, symbol="EP", mass=0, valent_radius=0.0,valents=[1]),# Ghost Atom
    Element(number=0, symbol="Bq", mass=0, valent_radius=0.0),  # Ghost Atom
    Element(number=0, symbol="D", mass=0, valent_radius=0.0),  # Ghost Atom
    Element(number=0, symbol="UNK", mass=0, valent_radius=0.0),  # Ghost Atom
    Element(number=0, symbol="DP", mass=0, valent_radius=0.0),  # Ghost Atom
    Element(number=0, symbol="VS", mass=0, valent_radius=0.0),  # Ghost Atom
]


_element_table = {}
for elem in ELEMENT_TABLE:
    _element_table[elem.symbol] = elem
    if elem.number > 0:
        _element_table[elem.number] = elem


def get_elem_property(label, target, value):
    if label == "mass":
        element = Element.get_by_mass(value)
    else:
        element = Element.get(value)

    if target in ["number","atom_number"]:
        return element.number
    elif target in ["element", "elem", "symbol"]:
        return element.symbol
    else:
        return getattr(element,target)

#def get_legacy_elements_table() -> Dict[str, list]:
#    return dict([element._to_legacy_tuple() for element in ALL_ELEMENTS])

def get_bonded_distance(e1, e2):
    elem1 = Element.get(e1)
    elem2 = Element.get(e2)
    return elem1.valent_radius + elem2.valent_radius

def get_bonded_type_distance(e1, e2, btype):
    elem1 = Element.get(e1)
    elem2 = Element.get(e2)
    v1 = elem1.valent_radius 
    v2 = elem2.valent_radius 
    if btype == "2":
        if elem1.double_bond_radius is not None:
            v1 = elem1.double_bond_radius
        else:
            return None
        if elem2.double_bond_radius is not None:
            v2 = elem2.double_bond_radius
        else:
            return None
    elif btype == "3":
        if elem1.triple_bond_radius is not None:
            v1 = elem1.triple_bond_radius
        else:
            return None
        if elem2.triple_bond_radius is not None:
            v2 = elem2.triple_bond_radius
        else:
            return None
    elif btype == "1.5" or btype == "ar":
        if elem1.conju_bond_radius is not None:
            v1 = elem1.conju_bond_radius
        else:
            return None
        if elem2.conju_bond_radius is not None:
            v2 = elem2.conju_bond_radius
        else:
            return None
    
    return v1 + v2

ELEMENT_ORDER = ["C","O","N","S","P","F","Cl","Br","I","B","Si","As","Te","Na","K","H"]
