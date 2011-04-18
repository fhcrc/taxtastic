from collections import defaultdict
from itertools import combinations

def union(it):
    ret = set()
    for x in it:
        ret.update(x)
    return ret

def walk(cur, metadata):
    parents, colors, cut_colors = metadata
    # The root node `T` has no `K(T)`, and choosing a `c` for that case makes
    # no sense. However, there still might be a `B`, so we fudge `K` to
    # preserve generality.
    if parents[cur] is None:
        K = None,
    else:
        K = cut_colors[cur]

    ret = defaultdict(dict)

    if not cur.clades:
        if K:
            color = colors[cur]
            assert len(K) == 1 and color in K
            ret[color][frozenset([color])] = {cur}
        else:
            ret[None] = {cur}
        return ret

    phi = [walk(x, metadata) for x in cur.clades]
    B = union(cut_colors[a] & cut_colors[b]
        for a, b in combinations(cur.clades, 2))

    # Trivial nodes are nodes below an edge with no cut colors.
    trivial_nodes = union(x.get(None, []) for x in phi)

    # If there were no cut colors, we don't need to do anything.
    if not K:
        ret[None] = trivial_nodes
        return ret

    # These phis will never matter, so there's no point in considering them.
    nontrivial_phis = [phi_i for phi_i in phi if not phi_i.get(None)]

    for c in K:
        ret_c = defaultdict(list)
        def aux(phis, used_colors, accum):
            # Base case; we've reached the end of the list.
            if not phis:
                ret_c[frozenset(used_colors)].append(accum)
                return

            phi_i, phi_rest = phis[0], phis[1:]
            for b in B | {c}:
                # The other part of the generality fudging mentioned above.
                # Since we always take the union of `B` and `{c}`, just ignore
                # the case where we pick the dummy `c`.
                if b is None:
                    continue
                X_is = phi_i[b]

                # One possible solution is to ignore this `phi` completely.
                aux(phi_rest, used_colors, accum)

                for X_i in phi_i[b]:
                    # Ignore a subcoloring if it's empty (which we've already
                    # handled) or it contains already-used colors (not counting
                    # `b`).
                    if not X_i or X_i & used_colors > {b}:
                        continue
                    aux(phi_rest, used_colors | X_i, accum | phi_i[b][X_i])

        aux(nontrivial_phis, set(), set())

        # For each `X_i`, the optimal `T_i` is the one with the most nodes in
        # it. Trivial nodes also will fit into any subtree, so add them back.
        ret[c] = {X_i: max(T_is, key=len) | trivial_nodes
            for X_i, T_is in ret_c.iteritems()}

    # Un-fudge the results if we're working on a root node.
    if parents[cur] is None:
        return ret[None]

    return ret
