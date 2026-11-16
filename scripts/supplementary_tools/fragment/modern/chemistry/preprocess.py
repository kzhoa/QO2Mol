from . import ring, structure


def prepare(m):
    elements, connections = structure.get_2d_connectivity(m)
    reduced, _, _ = structure.remove_hydrogen_atoms(elements, connections)
    cycles = ring.cyclo_property(m, ring.find_cyclo(reduced))
    rings, blocks, components = ring.cyclo_blocks(cycles)
    structure.assign_structure_info(m, rings, blocks, components)
    structure.assign_local_info(m, reduced)
    structure.assign_connect_type(m)
    structure.assign_conjugate_info(m)
    if len(m.Atoms) > 1:
        structure.divide_charge(m)
    m.steps.append("structure")
    structure.assign_hybrid(m)
    m.steps.append("hybrid")
    return m
