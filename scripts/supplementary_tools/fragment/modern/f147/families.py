"""147 family descriptors generated from an already partitioned graph."""

from copy import deepcopy
from .family_kernel import Fragmentation

FAMILIES = (
    "elem_frag",
    "seco_frag",
    "csf_frag",
    "rsf_frag",
    "chain_frag",
    "scaffold_frag",
    "tf_frag",
    "sssr_frag",
    "sketch_frag",
)

DEFAULT_FAMILIES = FAMILIES
BASELINE_ID = "modern_147_migrated_v1"


def describe(molecule, families=None):
    """Describe an already partitioned parent, preserving its atom indices."""
    selected = FAMILIES if families is None else tuple(families)
    unknown = set(selected) - set(FAMILIES)
    if unknown:
        raise ValueError("Unknown output families: " + str(sorted(unknown)))
    core = Fragmentation()
    core.run(molecule, selected)
    return {
        "profile": "legacy_147",
        "reconstruction_baseline": BASELINE_ID,
        "rule_count": 147,
        "elementary_schema": "elem_frag",
        "elementary_fragments": [
            {"index": i, "components": list(group)}
            for i, group in enumerate(molecule.frags)
        ],
        "fragments": {f: deepcopy(getattr(molecule, f)) for f in selected},
        "combination_status": "generated_from_parent_graph",
        "errors": [],
    }
