"""Modern product generation around the f147 kernel."""

from collections import Counter
from copy import deepcopy
from qqtools import LazyImport
from rdkit import rdBase
from .config import FragmentConfig
from .chemistry.graph import from_rdkit
from .chemistry.input import read_parent
from .chemistry.identity import identify
from .chemistry.descriptors import describe
from .chemistry.products import normalized_members, connected_members, cap


_BACKENDS = {"f147": LazyImport(f"{__package__}.f147")}
_f147_partition = LazyImport(f"{__package__}.f147.partition")
_f147_families = LazyImport(f"{__package__}.f147.families")


def validate_request(backend, families, include_parent, config):
    if backend != "f147":
        raise ValueError("Modern backend must be f147")
    if type(config) is not FragmentConfig:
        raise TypeError("config must be FragmentConfig")
    if type(include_parent) is not bool:
        raise TypeError("include_parent must be boolean")
    package = _BACKENDS[backend]
    if isinstance(families, str):
        raise TypeError("families must be a sequence of names, not a string")
    selected = package.DEFAULT_FAMILIES if families is None else tuple(families)
    if len(selected) != len(set(selected)) or set(selected) - set(package.FAMILIES):
        raise ValueError("Unknown or duplicate output families")
    return selected


def validate_partition(mol, fragments):
    rows = list(fragments.values()) if isinstance(fragments, dict) else fragments
    covered = []
    normalized = []
    for index, row in enumerate(rows):
        atoms = row["components"]
        if not atoms or any(
            type(i) is not int or not 0 <= i < mol.GetNumAtoms() for i in atoms
        ):
            raise ValueError("Invalid elementary atom membership")
        if len(connected_members(mol, atoms)) != 1:
            raise ValueError("Disconnected elementary partition")
        covered.extend(atoms)
        normalized.append(dict(index=index, components=sorted(atoms)))
    if Counter(covered) != Counter(range(mol.GetNumAtoms())):
        raise ValueError(
            "Elementary fragments must cover every parent atom exactly once"
        )
    return normalized


def process(record, families=None, include_parent=False, *, backend, config=None):
    config = FragmentConfig() if config is None else config
    selected = validate_request(backend, families, include_parent, config)
    mol, input_status = read_parent(record, config)
    parent = identify(mol)
    parent.update(
        identity_import=input_status, historical_inchi_key=record.get("inchi_key")
    )
    name = str(record.get("id", parent["identity_key"]))
    graph = from_rdkit(mol)
    graph = _f147_partition.partition(graph)
    seed = _f147_families.describe(graph, selected)
    adjustments = getattr(graph, "runtime_adjustments", [])
    elementary = validate_partition(mol, seed["elementary_fragments"])
    errors = [
        dict(e, stage="kernel_" + e.get("stage", "generation")) for e in seed["errors"]
    ]
    proposals = []
    for family, descriptors in seed["fragments"].items():
        for index, descriptor in descriptors.items():
            try:
                members = normalized_members(mol, descriptor["components"])
                for component_index, component in enumerate(
                    connected_members(mol, members)
                ):
                    proposals.append(
                        (family, index, component_index, component, descriptor)
                    )
            except Exception as exc:
                errors.append(
                    dict(
                        stage="membership",
                        family=family,
                        index=index,
                        error=f"{type(exc).__name__}: {exc}",
                    )
                )
            if len(proposals) + int(include_parent) > config.max_products_per_parent:
                raise ValueError(
                    "max_products_per_parent exceeded; no order-dependent truncation performed"
                )
    occurrences, candidates, excluded = [], {}, []

    def retain(identity, descriptors, provenance, construction):
        if identity["radical_electrons"] and not config.allow_radicals:
            raise ValueError("Unexpected radical product")
        heavy = descriptors["heavy_atoms"]
        if heavy < config.min_heavy_atoms or (
            config.max_heavy_atoms is not None and heavy > config.max_heavy_atoms
        ):
            excluded.append(
                dict(
                    provenance,
                    heavy_atoms=heavy,
                    identity_key=identity["identity_key"],
                    reason="outside_export_size_range",
                    cap_validation=construction["validation"],
                )
            )
            return
        key = identity["identity_key"]
        item = dict(
            identity,
            **provenance,
            construction=construction,
            descriptors=descriptors,
            classification_key=descriptors["classification_key"],
            classification_status="classified",
            parent_id=name,
            parent_identity_key=parent["identity_key"],
            source=record.get("source"),
            source_collection=record.get("collection"),
            source_scope=record.get("scope"),
        )
        if key not in candidates:
            candidates[key] = dict(
                identity, descriptors=descriptors, occurrence_indices=[], families=[]
            )
        candidates[key]["occurrence_indices"].append(len(occurrences))
        candidates[key]["families"] = sorted(
            set(candidates[key]["families"]) | {provenance["frag_type"]}
        )
        occurrences.append(item)

    # Cache only within this parent. Membership and internal bond edits fully
    # determine construction; per-occurrence provenance is added to fresh copies.
    product_cache = {}
    for family, index, component_index, members, descriptor in proposals:
        provenance = dict(
            frag_type=family,
            fragment_index=index,
            connected_component_index=component_index,
        )
        try:
            special = (
                tuple(sorted(set(descriptor.get("special_bond", ())) & set(members)))
                if family == "sssr_frag"
                else ()
            )
            key = (tuple(members), special)
            if key not in product_cache:
                product, construction = cap(mol, members, special_bond=special)
                product_cache[key] = (
                    identify(product),
                    describe(product),
                    construction,
                )
            identity, descriptors, construction = deepcopy(product_cache[key])
            construction["seed_components"] = list(descriptor["components"])
            construction["seed_offsets"] = descriptor.get("offsets", [])
            construction["retained_parent_hydrogens"] = sorted(
                set(members) - set(descriptor["components"])
            )
            construction["offsets_recomputed"] = sorted(construction["caps"]) != sorted(
                descriptor.get("offsets", [])
            )
            retain(identity, descriptors, provenance, construction)
        except Exception as exc:
            errors.append(
                dict(
                    provenance,
                    stage="product",
                    components=members,
                    error=f"{type(exc).__name__}: {exc}",
                )
            )
    if include_parent:
        retain(
            identify(mol),
            describe(mol),
            dict(
                frag_type="small_molecule",
                fragment_index=None,
                connected_component_index=None,
            ),
            dict(
                parent_atom_indices=list(range(mol.GetNumAtoms())),
                product_to_parent=list(range(mol.GetNumAtoms())),
                caps=[],
                hydrogen_caps=[],
                validation="validated_parent",
                normalization="rdkit_memory_boundary_v1",
                coordinates_generated=False,
            ),
        )
    return dict(
        id=name,
        backend=backend,
        profile=backend + "_v1",
        tool_version="1.2.0",
        parent=parent,
        elementary_fragments=elementary,
        seed_fragments=seed["fragments"],
        kernel_baseline=seed.get("reconstruction_baseline"),
        runtime_adjustments=adjustments,
        families=list(selected),
        config=config.to_dict(),
        config_sha256=config.fingerprint,
        rdkit_version=rdBase.rdkitVersion,
        occurrences=occurrences,
        candidates=[candidates[k] for k in sorted(candidates)],
        excluded=excluded,
        errors=errors,
        counts=dict(
            proposals=len(proposals) + int(include_parent),
            occurrences=len(occurrences),
            unique_candidates=len(candidates),
            excluded=len(excluded),
            errors=len(errors),
        ),
        selection="all_successful_exact_graph_identities",
        coordinates_generated=False,
        historical_reproduction=False,
        product_normalization="rdkit_memory_boundary_v1",
    )
