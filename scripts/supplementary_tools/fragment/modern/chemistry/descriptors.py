"""Complete classification calculated from the actual product graph."""

from collections import Counter
from rdkit import Chem
from rdkit.Chem import Lipinski


def describe(mol):
    heavy = Chem.RemoveHs(mol)
    rings = list(heavy.GetRingInfo().AtomRings())
    elements = dict(sorted(Counter(a.GetSymbol() for a in heavy.GetAtoms()).items()))
    charge = sum(a.GetFormalCharge() for a in heavy.GetAtoms())
    classification = dict(
        elements=sorted(elements),
        charge=charge,
        ring_sizes=sorted(len(r) for r in rings),
        aromatic_rings=sum(
            all(heavy.GetAtomWithIdx(i).GetIsAromatic() for i in r) for r in rings
        ),
    )
    return dict(
        heavy_atoms=sum(a.GetAtomicNum() > 1 for a in heavy.GetAtoms()),
        element_counts=elements,
        connected_components=len(Chem.GetMolFrags(heavy)),
        rotatable_bonds=Lipinski.NumRotatableBonds(heavy),
        classification_key=classification,
        classification_status="classified",
    )
