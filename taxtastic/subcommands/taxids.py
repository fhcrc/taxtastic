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
"""Convert a list of taxonomic names into a recursive list of species
leve tax_ids.

``The names to convert can be specified in a text file with one name
per line (the ``-f`` or ``--name-file`` options) or on the command
line as a comma delimited list (the ``-n`` of ``--name`` options).
"""
import logging
import argparse
import taxtastic
import sqlalchemy
import sys

from taxtastic.taxonomy import Taxonomy

log = logging.getLogger(__name__)


def get_children(engine, parent_ids, unordered_ranks, rank='species', schema=None):
    """
    Recursively fetch children of tax_ids in `parent_ids` until the
    rank of `rank`
    """

    if not parent_ids:
        return []

    nodes = schema + '.nodes' if schema else 'nodes'
    names = schema + '.names' if schema else 'names'

    cmd = ('select tax_id, tax_name, rank '
           'from {} join {} using (tax_id) '
           'where parent_id = :tax_id and is_primary').format(nodes, names)

    species = []
    for parent_id in parent_ids:
        result = engine.execute(sqlalchemy.sql.text(cmd), tax_id=parent_id)
        keys = list(result.keys())
        rows = [dict(list(zip(keys, row))) for row in result.fetchall()]
        for r in rows:
            if r['rank'] == rank and 'sp.' not in r['tax_name']:
                species.append(r)
        others = [r for r in rows
                  if r['rank'] != rank and r['rank'] not in unordered_ranks]
        if others:
            _, s = get_children(
                engine, [r['tax_id'] for r in others], unordered_ranks)
            species.extend(s)

    return keys, species


def build_parser(parser):
    parser = taxtastic.utils.add_database_args(parser)

    input_group = parser.add_argument_group(
        "Input options").add_mutually_exclusive_group()

    input_group.add_argument(
        '-f', '--name-file', metavar='FILE', type=argparse.FileType('rU'),
        dest='taxnames_file',
        help='file containing a list of taxonomic names, one per line')

    input_group.add_argument(
        '-n', '--name', metavar='NAMES',
        dest='taxnames',
        help=('list of taxonomic names provided as a comma-delimited '
              'list on the command line'))

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()
    output_group.add_argument(
        '-o', '--out', metavar='FILE', type=argparse.FileType('w'),
        help='output file', default=sys.stdout)


def action(args):

    taxnames_file = args.taxnames_file
    taxnames = args.taxnames

    engine = sqlalchemy.create_engine(args.url, echo=False)
    tax = Taxonomy(engine, schema=args.schema)

    names = []
    if taxnames_file:
        names += [line.split('#', 1)[0].strip()
                  for line in taxnames_file
                  if line.strip() and not line.startswith('#')]

    if taxnames:
        names += [x.strip() for x in taxnames.split(',')]

    taxa = {}
    for name in set(names):
        tax_id, tax_name, is_primary, rank, note = '', '', '', '', ''

        try:
            tax_id, tax_name, is_primary = tax.primary_from_name(name)
        except ValueError:
            log.warning(name + ' not found')
            continue
        else:
            parent, rank = tax._node(tax_id)
            note = '' if is_primary else 'not primary'

        if note:
            log.warning(
                '%(name)20s | %(tax_id)7s %(tax_name)20s %(note)s' % locals())

        if rank == 'species':
            taxa[tax_id] = dict(tax_id=tax_id, tax_name=tax_name, rank=rank)
        else:
            keys, rows = get_children(
                engine, [tax_id], tax.unordered_ranks, schema=args.schema)
            taxa.update(dict((row['tax_id'], row) for row in rows))

    for d in sorted(list(taxa.values()), key=lambda x: x['tax_name']):
        args.out.write('%(tax_id)s # %(tax_name)s\n' % d)
