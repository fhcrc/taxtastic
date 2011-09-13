# This file is part of taxtastic.
#
#    taxtastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    taxtastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.
import collections
import logging
from itertools import combinations

log = logging.getLogger(__name__)

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

    # The root node is reported to cut every color that crosses the root, but
    # we want to treat the root node as if it doesn't cut any color because we
    # only want the best set of nodes from the root.
    if parents[cur] is None:
        K = set()
    else:
        K = cut_colors[cur]

    ret = collections.defaultdict(dict)

    if not cur.clades:
        if K:
            color = colors[cur]
            assert len(K) == 1 and color in K
            ret[color][frozenset([color])] = {cur}
            ret[None][frozenset([color])] = {cur}
        else:
            ret[None][frozenset()] = {cur}
        return ret

    phi = [walk(x, metadata) for x in cur.clades]
    B = union(cut_colors[a] & cut_colors[b]
        for a, b in combinations(cur.clades, 2))

    for c in K | {None}:
        ret_c = collections.defaultdict(list)
        for b in B | {c}:
            def aux(phis, used_colors, accum):
                # Base case; we've reached the end of the list.
                if not phis:
                    ret_c[frozenset(used_colors)].append(accum)
                    return

                phi_i, phi_rest = phis[0], phis[1:]
                X_is = phi_i[b]

                # One possible solution is to ignore this `phi` completely.
                aux(phi_rest, used_colors, accum)

                if not X_is:
                    X_is = phi_i[None]

                for X_i in X_is:
                    if (X_i & used_colors) - {b}:
                        continue
                    if b != c and c in X_i:
                        continue
                    aux(phi_rest, used_colors | X_i, accum | X_is[X_i])

            aux(phi, set(), set())

        # For each `X_i`, the optimal `T_i` is the one with the most nodes in
        # it.
        ret[c] = {X_i: max(T_is, key=len) for X_i, T_is in ret_c.iteritems()}

    # If there were no cut colors, the only relevant data is the biggest set of
    # nodes, so prune everything else out.
    if not K:
        total = max(ret[None].itervalues(), key=len)
        ret.clear()
        ret[None][frozenset()] = total

    # If this is the parent node, return just the biggest set of nodes.
    if parents[cur] is None:
        return ret[None][frozenset()]

    return ret

Ranking = collections.namedtuple('Ranking', 'rank node')

def reroot_from_rp(root, rp, ignore_missing_sequences=False):
    name_map = dict(rp.db.cursor().execute("""
        SELECT seqname, tax_id
        FROM   sequences
    """))
    rank_map = dict(rp.db.cursor().execute("""
        SELECT tax_id, rank_order
        FROM   taxa
               JOIN ranks USING (rank)
    """))
    def subrk_min(t):
        mrca = rp.most_recent_common_ancestor(
            *set(name_map[n.name] for n in t.get_terminals()
                 if not ignore_missing_sequences or n.name in name_map))
        logging.debug("mrca for %r is %r", t, mrca)
        return rank_map[mrca]

    return reroot(root, subrk_min)

def reroot(cur, subrk_min):
    while True:
        if len(cur.clades) < 2:
            return cur

        ranks = sorted(Ranking(subrk_min(n), n) for n in cur)
        logging.debug("rankings: %r", ranks)
        if ranks[0].rank == ranks[1].rank:
            return cur
        cur = ranks[0].node
