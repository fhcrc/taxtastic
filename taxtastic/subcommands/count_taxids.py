"""Count tax_id appearances in a taxtable lineage

Can narrow tax_id count by seq_info if provided.

Returns taxonomy with columns ['tax_id', 'tax_name', 'rank', 'count']
"""
import pandas
import sys


def build_parser(p):
    # inputs
    p.add_argument(
        'taxonomy',
        help="""Taxonomy metadata, minimum columns with
        ordered column taxonomy - [tax_id, root,...,{last_rank}]""")
    p.add_argument(
        '-i', '--seq-info',
        help=('csv of actual sequence representatives. '
              'Minimum column [tax_id]'))
    p.add_argument(
        '--out',
        default=sys.stdout,
        help=('taxonomy output with counts'))


def action(args):
    lineage = pandas.read_csv(args.taxonomy, dtype=str)

    if args.seq_info:
        seqinfo = pandas.read_csv(args.seq_info, usecols=['tax_id'], dtype=str)
        taxonomy = seqinfo.merge(lineage, on='tax_id', how='left')
    else:
        taxonomy = lineage

    tax_cols = taxonomy.columns.tolist()
    ranks = tax_cols[tax_cols.index('root'):]
    rank_counts = [taxonomy[r].value_counts(sort=False) for r in ranks]

    counts = pandas.concat(rank_counts)
    counts.name = 'count'

    results = lineage.join(counts, on='tax_id', how='right')

    # apply rank indexing to sort by rank
    def rank_index(row):
        return ranks.index(row['rank'])
    results['rank_index'] = results.apply(rank_index, axis=1)

    # sort by [count, rank, tax_name] in that priority
    results = results.sort_values(
        by=['count', 'rank_index', 'tax_name'], ascending=[False, True, True])

    out_cols = ['tax_id', 'tax_name', 'rank', 'count']
    results[out_cols].to_csv(args.out, index=False)
