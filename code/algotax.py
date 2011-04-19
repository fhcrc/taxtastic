import collections
from itertools import combinations

def union(it):
    "Take the union of an iterable of sets."

    ret = set()
    for x in it:
        ret.update(x)
    return ret

def intersection(it):
    "Take the intersection of an iterable of sets."

    ret = None
    for x in it:
        if ret is None:
            ret = set(x)
        else:
            ret.intersection_update(x)
    return ret or set()

CladeMetadata = collections.namedtuple(
    'CladeMetadata', 'parents colors cut_colors')

def color_clades(tree, colors):
    "Given a biopython tree and colors of its leaves, color its edges."

    parents = {tree.root: None}
    cut_colors = collections.defaultdict(set)
    stack = [('down', tree.root, None)]
    while stack:
        phase, cur, color = stack.pop()
        if phase == 'down':
            if not cur.clades:
                if cur not in colors:
                    continue
                stack.append(('up', cur, colors[cur]))
            else:
                for child in cur.clades:
                    parents[child] = cur
                    stack.append(('down', child, None))
        elif phase == 'up':
            if cur is None or color in cut_colors[cur]:
                continue
            cut_colors[cur].add(color)
            stack.append(('up', parents[cur], color))

    stack = [(tree.root, set())]
    while stack:
        node, okayed = stack.pop()
        if not node.clades:
            continue
        okayed = union(cut_colors[a] & cut_colors[b]
            for a, b in combinations(node.clades, 2)) | okayed
        for e in node.clades:
            e_ = cut_colors[e] & okayed
            if e_ != cut_colors[e]:
                stack.append((e, okayed))
                cut_colors[e] = e_
    return CladeMetadata(parents, colors, cut_colors)

def walk(cur, metadata):
    "Walk a biopython clade, determining the optimal convex subcoloring."

    parents, colors, cut_colors = metadata
    # Choosing None for `c` basically means 'coalesce disjoint colorings',
    # which is what we want to do at the root node.
    if parents[cur] is None:
        K = None,
    else:
        K = cut_colors[cur]

    ret = collections.defaultdict(dict)

    if not cur.clades:
        if K:
            color = colors[cur]
            assert len(K) == 1 and color in K
            ret[color][frozenset([color])] = {cur}
        else:
            ret[None][frozenset()] = {cur}
        return ret

    phi = [walk(x, metadata) for x in cur.clades]
    B = union(cut_colors[a] & cut_colors[b]
        for a, b in combinations(cur.clades, 2))

    # As above, we want to coalesce the colorings. Since this edge isn't
    # colored, the different possibilities don't matter to edges above it; just
    # choose the best possibility.
    if not K:
        K = None,

    for c in K:
        ret_c = collections.defaultdict(list)
        def aux(phis, used_colors, accum):
            # Base case; we've reached the end of the list.
            if not phis:
                ret_c[frozenset(used_colors)].append(accum)
                return

            phi_i, phi_rest = phis[0], phis[1:]
            for b in B | {c, None}:
                X_is = phi_i[b]

                # One possible solution is to ignore this `phi` completely.
                aux(phi_rest, used_colors, accum)

                for X_i in phi_i[b]:
                    # Ignore a subcoloring if it's empty (which we've already
                    # handled) or it contains already-used colors (not counting
                    # `b`).
                    if b is not None and X_i & used_colors > {b}:
                        continue
                    aux(phi_rest, used_colors | X_i, accum | phi_i[b][X_i])

        aux(phi, set(), set())

        # For each `X_i`, the optimal `T_i` is the one with the most nodes in
        # it.
        ret[c] = {X_i: max(T_is, key=len) for X_i, T_is in ret_c.iteritems()}

    # If this is the parent node, it's more useful to be back the coalesced
    # results than something mapping to them.
    if parents[cur] is None:
        return ret[None]

    # Otherwise if there were no cut colors, the only relevant data is the
    # biggest set of nodes, so prune everything else out.
    elif ret.get(None):
        total = max(ret[None].itervalues(), key=len)
        ret.clear()
        ret[None][frozenset()] = total

    return ret
