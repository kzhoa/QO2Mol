"""Graph/structure helpers shared with the sibling fragment toolkit."""

if "." in __package__:
    from .. import fragment as fragment_package
    from ..fragment.common.elements import Element
    from ..fragment.common.graph import Molecule, from_document, from_rdkit
    from ..fragment.common.preprocess import prepare
else:
    import fragment as fragment_package
    from fragment.common.elements import Element
    from fragment.common.graph import Molecule, from_document, from_rdkit
    from fragment.common.preprocess import prepare
