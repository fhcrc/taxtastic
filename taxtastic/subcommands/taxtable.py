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
import pandas

import re

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import getlines

from sqlalchemy import create_engine
import os.path
import sys

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'database_file',
        metavar='SQLITE',
        help='Name of the sqlite database file')

    parser.add_argument(
        '--from-table',
        metavar='CSV',
        help=('taxtable to derive new taxtable from'))
    parser.add_argument(
        '--full',
        action='store_true',
        help='Include rank columns not in final lineages.')

    input_group = parser.add_argument_group('input options')

    input_group.add_argument(
        '-n', '--tax-names',
        dest='taxnames',
        metavar='FILE',
        help="""A file identifing taxa in the form of taxonomic
        names. Names are matched against both primary names and
        synonyms. Lines beginning with "#" are ignored. Taxa
        identified here will be added to those specified using
        --tax-ids""")

    input_group.add_argument(
        '-t', '--tax-ids',
        dest='taxids',
        metavar='FILE-OR-LIST',
        help="""File containing a whitespace-delimited list of
        tax_ids (ie, separated by tabs, spaces, or newlines; lines
        beginning with "#" are ignored). This option can also be
        passed a comma-delited list of taxids on the command line.""")

    input_group.add_argument(
        '-i', '--seq-info', type=argparse.FileType('r'),
        help="""Read tax_ids from sequence info file, minimally containing the
        field "tax_id" """)

    input_group.add_argument(
        '--from-id',
        help=('taxid to branch from'))

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-o', '--out-file',
        dest='out_file',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='FILE',
        help="""Output file containing lineages for the specified taxa
        in csv format; writes to stdout if unspecified""")


def action(args):
    engine = create_engine(
        'sqlite:///%s' % args.database_file, echo=args.verbosity > 2)

    ranks = pandas.read_sql_table('ranks', engine, index_col='index')
    ranks = ranks['rank'].tolist()

    if args.from_table:
        log.info('using taxtable ' + args.from_table)
        taxtable = pandas.read_csv(args.from_table, dtype=str).set_index('tax_id')
    else:
        log.info('building taxtable from ' + args.database_file)
        nodes = pandas.read_sql_table('nodes', engine, index_col='tax_id')
        names = pandas.read_sql_table(
            'names', engine, columns=['tax_id', 'tax_name', 'is_primary'])
        names = names[names['is_primary']].set_index('tax_id')
        nodes = pandas.read_sql_table('nodes', engine, index_col='tax_id')
        nodes = nodes.join(names['tax_name'])
        taxtable = build_taxtable(nodes, ranks)

    if args.from_id:
        from_taxon = taxtable.loc[args.from_id]

        # select all rows where rank column == args.from_id
        from_table = taxtable[taxtable[from_taxon['rank']] == args.from_id]

        # build taxtable up to root from args.from_id
        while from_taxon.name != '1':  # root
            parent = taxtable.loc[from_taxon['parent_id']]
            from_table = pandas.concat([pandas.DataFrame(parent).T, from_table])
            from_taxon = parent
        # reset lost index name after concatenating transposed series
        from_table.index.name = 'tax_id'
        taxtable = from_table

    if any([args.taxids, args.taxnames, args.seq_info]):
        tax = Taxonomy(engine, ranks)
        tax_ids = set()
        if args.taxids:
            if os.access(args.taxids, os.F_OK):
                for line in getlines(args.taxids):
                    tax_ids.update(set(re.split(r'[\s,;]+', line)))
            else:
                tax_ids.update(
                    [x.strip() for x in re.split(r'[\s,;]+', args.taxids)])

        if args.seq_info:
            with args.seq_info:
                reader = csv.DictReader(args.seq_info)
                tax_ids.update(
                    frozenset(i['tax_id'] for i in reader if i['tax_id']))

        if not(are_valid(tax_ids, tax)):
            return "Some taxids were invalid.  Exiting."

        if args.taxnames:
            for taxname in getlines(args.taxnames):
                for name in re.split(r'\s*[,;]\s*', taxname):
                    tax_id, primary_name, is_primary = tax.primary_from_name(
                        name.strip())
                    tax_ids.add(tax_id)

        keepers = taxtable.loc[tax_ids]
        for col in keepers.columns:
            if col in ranks:
                tax_ids.update(keepers[col].dropna().values)
        taxtable = taxtable.loc[tax_ids]

    if not args.full:
        taxtable = taxtable.dropna(axis=1, how='all')

    # sort columns
    taxtable = taxtable[
        ['rank', 'tax_name'] + [r for r in ranks if r in taxtable.columns]]

    # sort rows
    taxtable['rank'] = taxtable['rank'].astype('category', categories=ranks)
    taxtable = taxtable.sort_values('rank')

    # write and close db
    taxtable.to_csv(args.out_file)
    engine.dispose()


def are_valid(tax_ids, tax):
    '''
    Check if ALL tax_ids are valid.  Return True/False
    '''
    valid = True
    for t in tax_ids:
        try:
            tax._node(t)
        except ValueError:
            # Check for merged
            m = tax._get_merged(t)
            if m and m != t:
                msg = ("Taxid {0} has been replaced by {1}. "
                       "Please update your records").format(t, m)
                print >> sys.stderr, msg
            else:
                print >>sys.stderr, "Taxid %s not found in taxonomy." % t
            valid = False

    return valid


def build_taxtable(df, ranks):
    '''
    given list of tax_ids with parent_ids and an ordered list of ranks return
    a table of taxonomic lineages with ranks as columns
    '''

    df_index_name = df.index.name
    rank_count = len(ranks)

    df = df.join(df['rank'], on='parent_id', rsuffix='_parent').reset_index()
    df['rank_parent'] = df['rank_parent'].astype('category', categories=ranks)
    lineages = df[df['tax_id'] == '1'].iloc[[0]]
    lineages.loc[:, 'root'] = lineages['tax_id']
    df = df.drop(lineages.index)
    tax_parent_groups = df.groupby(by='rank_parent', sort=True)

    for i, (parent, pdf) in enumerate(tax_parent_groups):
        at_rank = []
        for child, df in pdf.groupby(by='rank', sort=False):
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

        sys.stderr.write('{} of {} rank lineages completed\r'.format(i, rank_count))

    return lineages.drop('rank_parent', axis=1).set_index(df_index_name)
