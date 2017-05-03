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
"""Create a tabular representation of taxonomic lineages

Write a CSV file containing the minimal subset of the taxonomy in
``database_file`` which encompasses all the taxa specified in
``taxa_names.txt`` and ``tax_ids`` and all nodes connecting them to
the root of the taxonomy.  By default the CSV is written to
``stdout``, unless redirectored with ``-o`` or ``--out-file``.

``taxa_names.txt`` should be a text file specifying names of taxa.
Python style comments are ignored as are empty lines.  Names may be
separated by commas, semicolons, and arbitrary amounts of whitespace
on both sides of those separators, but the whitespace within them must
be exact (e.g., ``Lactobacillus crispatus`` must have exactly one space
between the two words to match the entry in the taxonomy).

``tax_ids`` is either a comma or semicolon delimited list of tax_ids
(e.g., ``4522,2213;44;221``) or the name of a text file containing
tax_ids.  The text file also allows Python style comments, and any
non-comment text separated by and combination of spaces, commas, and
semicolons is considered a tax_id.

tax_ids and taxa names can overlap, nor does anything have to be
unique within either file.  The nodes will only be written once in the
CSV output no matter how many times a particular taxon is mentioned.
"""
import argparse
import csv
import logging
import os
import pandas
import re
import sqlalchemy
import sys

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import getlines, add_database_args

log = logging.getLogger(__name__)


def build_parser(parser):
    parser = add_database_args(parser)

    node_parser = parser.add_argument_group(title='node options')
    node_parser.add_argument(
        '--valid',
        action='store_true',
        help='include only valid nodes')
    node_parser.add_argument(
        '--ranked',
        action='store_true',
        help='include only ranked nodes (drop no_rank nodes)')

    input_group = parser.add_argument_group('input options')

    input_group.add_argument(
        '--taxtable',
        metavar='CSV',
        help='build from a previous built taxtable')

    input_group.add_argument(
        '--clade-ids',
        help=('return top-down tax_id clades'))

    input_group.add_argument(
        '-n', '--taxnames',
        metavar='FILE',
        help=('A file identifing taxa in the form of taxonomic '
              'names. Names are matched against both primary names and '
              'synonyms. Lines beginning with "#" are ignored. Taxa '
              'identified here will be added to those specified using '
              '--tax-ids'))

    input_group.add_argument(
        '-t', '--tax-ids',
        metavar='FILE-OR-LIST',
        help=('File containing a whitespace-delimited list of '
              'tax_ids (ie, separated by tabs, spaces, or newlines; lines '
              'beginning with "#" are ignored). This option can also be '
              'passed a comma-delited list of tax_ids on the command line.'))

    input_group.add_argument(
        '-i', '--seq-info',
        type=argparse.FileType('r'),
        help=('Read tax_ids from sequence info file, minimally '
              'containing the field "tax_id"'))

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-o', '--out',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='FILE',
        help=('Output file containing lineages for the specified taxa '
              'in csv format; writes to stdout if unspecified'))


def action(args):
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)

    ranks_df = pandas.read_sql_table(
        'ranks', engine,
        schema=args.schema)
    # most operations in this script require ordering from 'root' down
    ranks_df = ranks_df.sort_values(by='height', ascending=False)
    ranks = ranks_df['rank'].tolist()

    nodes = None
    subset_ids = set()

    # check tax_ids subsets first before building taxtable
    if any([args.tax_ids, args.taxnames, args.seq_info]):
        tax = Taxonomy(engine, schema=args.schema)
        if args.tax_ids:
            if os.access(args.tax_ids, os.F_OK):
                for line in getlines(args.tax_ids):
                    subset_ids.update(set(re.split(r'[\s,;]+', line)))
            else:
                subset_ids.update(
                    [x.strip() for x in re.split(r'[\s,;]+', args.tax_ids)])

        if args.seq_info:
            log.info('reading tax_ids ' + args.seq_info.name)
            with args.seq_info:
                reader = csv.DictReader(args.seq_info)
                subset_ids.update(
                    frozenset(i['tax_id'] for i in reader if i['tax_id']))

        # this will raise an error if any tax_ids do not exist in database
        all_known(subset_ids, tax)

        if args.taxnames:
            for taxname in getlines(args.taxnames):
                for name in re.split(r'\s*[,;]\s*', taxname):
                    tax_id, primary_name, is_primary = tax.primary_from_name(
                        name.strip())
                    subset_ids.add(tax_id)

        if not subset_ids:
            log.error('no tax_ids to subset taxtable, exiting')
            return

    log.info('loading nodes table from database')
    nodes = pandas.read_sql_table(
        'nodes', engine, schema=args.schema, index_col='tax_id')

    if args.taxtable:
        log.info('using existing taxtable ' + args.taxtable)
        taxtable = pandas.read_csv(args.taxtable, dtype=str)
        taxtable = taxtable.set_index('tax_id')
        taxtable = taxtable.join(nodes[['parent_id', 'is_valid']])
    else:
        log.info('building taxtable')
        names = pandas.read_sql_table(
            'names', engine,
            schema=args.schema,
            columns=['tax_id', 'tax_name', 'is_primary'])
        names = names[names['is_primary']].set_index('tax_id')
        len_nodes = len(nodes)
        nodes = nodes.join(names['tax_name'])
        assert len_nodes == len(nodes)
        taxtable = build_taxtable(nodes, ranks)

    # subset taxtable clade lineages
    if args.clade_ids:
        dtypes = taxtable.dtypes
        clades = []
        for i in args.clade_ids.split(','):
            ancestor = taxtable.loc[i]

            # select all rows where rank column == args.from_id
            clade = taxtable[taxtable[ancestor['rank']] == i]

            # build taxtable up to root from args.from_id
            while ancestor.name != '1':  # root
                parent = taxtable.loc[ancestor['parent_id']]
                clade = pandas.concat(
                    [pandas.DataFrame(parent).T, clade])
                ancestor = parent
            # reset lost index name after concatenating transposed series
            clades.append(clade)
        taxtable = pandas.concat(clades)
        taxtable = taxtable[~taxtable.index.duplicated()]

        # set index.name and dtypes back after concating transposed series
        taxtable.index.name = 'tax_id'
        for d, t in dtypes.iteritems():
            taxtable[d] = taxtable[d].astype(t)

    # subset taxtable by set of tax_ids
    if subset_ids:
        keepers = taxtable.loc[subset_ids]
        for col in keepers.columns:
            if col in ranks:
                subset_ids.update(keepers[col].dropna().values)
        taxtable = taxtable.loc[subset_ids]

    # drop no rank nodes
    if args.ranked:
        ranks = ranks_df[~ranks_df['no_rank']]['rank'].tolist()
        taxtable = taxtable[taxtable['rank'].isin(ranks)]

    if args.valid:
        invalid = taxtable[~taxtable['is_valid']]
        # remove all invalids from the rank columns
        for r, g in invalid.groupby(by='rank'):
            taxtable.loc[taxtable[r].isin(g.index), r] = None
        # remove invalid rows
        taxtable = taxtable[taxtable['is_valid']]

    # clean up empty rank columns
    taxtable = taxtable.dropna(axis=1, how='all')

    # sort final column output
    taxtable = taxtable[
        ['rank', 'tax_name'] + [r for r in ranks if r in taxtable.columns]]

    # sort rows
    taxtable['rank'] = taxtable['rank'].astype('category', categories=ranks)
    taxtable = taxtable.sort_values('rank')

    # write and close db
    taxtable.to_csv(args.out)
    engine.dispose()


def all_known(tax_ids, tax):
    '''
    Check if ALL tax_ids are known.  Return True/False
    '''
    all_known = True
    for t in tax_ids:
        try:
            tax._node(t)
        except ValueError:
            # Check for merged
            m = tax._get_merged(t)
            if m and m != t:
                msg = ("Taxid {} has been replaced by {}. "
                       "Please update your records").format(t, m)
                log.error(msg)
            else:
                log.error('Taxid {} not found in taxonomy'.format(t))
            all_known = False

    if not all_known:
        raise ValueError('Some tax_ids are unknown.  Exiting.')


def build_taxtable(nodes, ranks):
    '''
    Given list of tax_ids with parent_ids and an ordered list of ranks return
    a table of taxonomic lineages with ranks as columns

    Iterate by rank starting with root and joining on parent_id to build
    lineages.
    '''
    df_index_name = nodes.index.name

    nodes = nodes.join(nodes['rank'], on='parent_id', rsuffix='_parent').reset_index()
    nodes['rank_parent'] = nodes['rank_parent'].astype('category', categories=ranks)
    lineages = nodes[nodes['tax_id'] == '1'].iloc[[0]]
    lineages.loc[:, 'root'] = lineages['tax_id']
    nodes = nodes.drop(lineages.index)
    tax_parent_groups = nodes.groupby(by='rank_parent', sort=True)
    # standard message width so everything prints nice
    msg = 'Processing lineages: {:<' + str(max(map(len, ranks))) + '}\r'
    for parent, pdf in tax_parent_groups:
        at_rank = []
        for child, df in pdf.groupby(by='rank'):
            df = df.copy()
            df[child] = df['tax_id']
            at_rank.append(df)
        if at_rank:
            '''bug in pandas 0.18.1 puts all categories as group
               keys even when there are no values so we much test
               for at_rank to be not empty'''
            at_rank = pandas.concat(at_rank)
            at_rank = at_rank.merge(
                lineages[ranks[:ranks.index(parent) + 1]],
                left_on='parent_id',
                right_on=parent,
                how='inner')
            lineages = lineages.append(at_rank)

        # status message
        sys.stderr.write(msg.format(parent))

    return lineages.drop('rank_parent', axis=1).set_index(df_index_name)
