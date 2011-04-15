def union(it):
    ret = set()
    for x in it:
        ret.update(x)
    return ret

def walk(cur):
    K = cur.root_edge.colors
    if cur.is_leaf:
        ret = defaultdict(dict)
        for color in K:
            ret[color][{cur.color}] = {cur}
            ret[color][{}] = {}
        if not ret:
            ret[None] = {cur}
        return ret

    phi = [walk(x.bottom) for x in cur.children]
    B = union(a.colors & b.colors
        for a, b in itertools.combinations(cur.children, 2))

    ret = defaultdict(dict)
    trivial_nodes = union(x[None] for x in phi)
    for c in K:
        ret_c = ret[c]
        def aux(phis, used_colors, accum):
            if not phis:
                ret_c[frozenset(used_colors)] = accum
            phi_i, phi_rest = phis[0], phis[1:]
            ret = []
            for b in B | {c}:
                for X_i in phi_i[b]:
                    if X_i & used_colors > {c}:
                        continue
                    aux(phi_rest, used_colors | X_i, accum | phi_i[b][X_i])
        aux(phi, set(), set())
    if not ret:
        ret[None] = trivial_nodes
    return ret
