"""Validated elementary partition -> real boundary caps -> strict identities."""
from dataclasses import asdict
from rdkit import rdBase
from .models import QQConfig, VERSION, TOOL_VERSION, FAMILIES, implementation_hash
from .input import read_parent
from .boundaries import partition_boundaries
from .products import cap
from .identity import identify
from .descriptors import describe
from .modules import describe_module


def validate_request(families, include_parent, config):
    if type(config) is not QQConfig:
        raise TypeError('config must be QQConfig from this release')
    if type(include_parent) is not bool or include_parent:
        raise ValueError('include_parent=True is removed; production outputs elementary only. See README.md')
    selected = FAMILIES if families is None else tuple(families)
    if selected != FAMILIES:
        raise ValueError('Only families=("elementary",) is supported; combinations/context are not supported. See README.md')


def process(record, families=None, include_parent=False, config=None):
    config = QQConfig() if config is None else config
    validate_request(families, include_parent, config)
    mol, input_status = read_parent(record, config)
    parent = identify(mol)
    name = str(record.get('id', parent['identity_key']))
    parent.update(identity_import=input_status, historical_inchi_key=record.get('inchi_key'))
    parent['historical_identity_matches'] = (None if not record.get('inchi_key') else record['inchi_key'] == parent['inchi_key'])
    (elementary, adjacency, decisions, motifs, operations), boundary_operations = partition_boundaries(mol, config)
    allowed_ring_bonds = [mol.GetBondBetweenAtoms(*edge).GetIdx()
                          for edge in boundary_operations["macrocycles"]["new_ring_cuts"]]
    if len(elementary) > config.max_candidates:
        raise ValueError('QQ candidate limit exceeded; no order-dependent truncation performed')
    occurrences, excluded, errors, elementary_output = [], [], [], []
    candidates = {}
    for fragment in elementary:
        item = asdict(fragment)
        heavy = len(fragment.heavy_atoms)
        item['within_candidate_size_range'] = config.min_heavy_atoms <= heavy <= config.max_heavy_atoms
        proposal = dict(fragment_indices=[fragment.index], components=list(fragment.atoms), origins=[{'family':'elementary'}])
        if not item['within_candidate_size_range']:
            excluded.append(dict(proposal, reason='too_small' if heavy < config.min_heavy_atoms else 'too_large', heavy_atoms=heavy))
        # Every block is capped/audited even outside export limits; export
        # filtering cannot hide large elementary blocks or their errors.
        try:
            product, construction = cap(mol, fragment.atoms, allowed_ring_bonds=allowed_ring_bonds)
            identity = identify(product)
            assembly = describe_module(product, construction)
            item.update(cap_status='ok', assembly=assembly, construction=construction,
                        capped_identity_key=identity['identity_key'])
            if item['within_candidate_size_range']:
                key = identity['identity_key']
                occurrence = dict(candidate_key=key, parent_id=name, parent_identity_key=parent['identity_key'],
                    source=record.get('source'), origins=proposal['origins'], fragment_indices=[fragment.index],
                    construction=construction, assembly=assembly)
                if key not in candidates:
                    candidates[key] = dict(identity, descriptors=describe(product, config),
                        occurrence_indices=[], families=['elementary'], assembly_keys=[])
                index = len(occurrences)
                occurrences.append(occurrence)
                candidates[key]['occurrence_indices'].append(index)
                candidates[key]['assembly_keys'] = sorted(set(candidates[key]['assembly_keys']) | {assembly['module_key']})
        except Exception as exc:
            message = f'{type(exc).__name__}: {exc}'
            item.update(cap_status='failed', assembly_error=message)
            errors.append(dict(stage='elementary_product', fragment_index=fragment.index, components=list(fragment.atoms), error=message))
        elementary_output.append(item)
    return dict(id=name, backend='qq', profile=VERSION, tool_version=TOOL_VERSION,
        config=config.to_dict(), config_sha256=config.fingerprint, implementation_sha256=implementation_hash(),
        rdkit_version=rdBase.rdkitVersion, parent=parent, input_atom_count=mol.GetNumAtoms(),
        bond_decisions=[asdict(d) for d in decisions], structural_matches=motifs,
        basic_unit_operations=operations, boundary_operations=boundary_operations,
        elementary_fragments=elementary_output,
        fragment_adjacency={str(i):sorted(v) for i,v in adjacency.items()},
        occurrences=occurrences, candidates=[candidates[k] for k in sorted(candidates)], excluded=excluded, errors=errors,
        counts=dict(proposals=len(elementary), occurrences=len(occurrences), unique_candidates=len(candidates),
                    excluded=len(excluded), errors=len(errors)),
        selection='all_successful_strict_identities', coordinates_generated=False)
