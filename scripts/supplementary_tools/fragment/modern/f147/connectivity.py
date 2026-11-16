"""Connected-component helper for 147 grouping."""


def find_mole(connect_dict):
    """Return disjoint, sorted components without modifying the input graph.

    Only keyed vertices belong to the induced graph. Connections are undirected,
    including when the caller supplies only one direction of an edge.
    """
    adjacency = {a: set() for a in connect_dict}
    for a, neighbors in connect_dict.items():
        for b in neighbors:
            if b in adjacency:
                adjacency[a].add(b)
                adjacency[b].add(a)
    unseen = set(adjacency)
    components = []
    for start in sorted(adjacency):
        if start not in unseen:
            continue
        unseen.remove(start)
        pending, component = [start], []
        while pending:
            a = pending.pop()
            component.append(a)
            for b in sorted(adjacency[a]):
                if b in unseen:
                    unseen.remove(b)
                    pending.append(b)
        components.append(sorted(component))
    return components
