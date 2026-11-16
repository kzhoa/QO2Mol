"""Ordered local refinement of elementary molecular boundaries."""

from . import (
    conflicts,
    acyl_periphery,
    aryl_ether_periphery,
    aryl_amine_periphery,
    macrocycles,
)


def partition_boundaries(mol, config):
    partition, conflict_records = conflicts.resolve(mol, config)
    records = {"conflicts": conflict_records}
    for name, module in (
        ("acyl_periphery", acyl_periphery),
        ("aryl_ether_periphery", aryl_ether_periphery),
        ("aryl_amine_periphery", aryl_amine_periphery),
        ("macrocycles", macrocycles),
    ):
        partition, records[name] = module.resolve(mol, config, partition)
    return partition, records
