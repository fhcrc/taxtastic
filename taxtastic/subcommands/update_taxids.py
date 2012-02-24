#!/usr/bin/env python
"""
Update obsolete tax_ids in preparation for use in taxtable. Takes sequence info
file as passed to `taxit create --seq-info`
"""
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
import argparse
import csv
import logging
import os.path
import sys

import sqlalchemy

from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile', help="""Input CSV file to process, minimally
            containing the fields 'seqname' and 'tax_id'. Rows with missing
            tax_ids are left unchanged.""", type=argparse.FileType('r'))
    parser.add_argument('-d', '--database-file', help="""Path to the taxonomy
            database""", required=True)
    parser.add_argument('-o', '--out-file', help="""Output file to write
            updates [default: stdout]""", type=argparse.FileType('w'),
            default=sys.stdout)
    parser.add_argument('-u', '--unknown-action', help="""Action to take on
            encountering an unknown tax_id [default: %(default)s]""",
            choices=('halt', 'remove'), default='halt')


def load_csv(fp):
    """
    Load a CSV file from fp, returning (headers, dialect, row_iter)

    fp must be seekable.
    """
    sample = fp.read(1024)
    fp.seek(0)
    dialect = csv.Sniffer().sniff(sample)
    reader = csv.DictReader(fp, dialect=dialect)
    return (reader.fieldnames, dialect, reader)

def taxid_updater(taxonomy, action='halt'):
    """
    Create a function to update obsolete taxonomies in
    """
    def update_taxid(row):
        current_tax_id = row['tax_id']
        if not current_tax_id:
            # Skip blank
            return row
        try:
            node = taxonomy._node(current_tax_id)
            # _node raises KeyError if the taxon couldn't be found in the
            # current taxonomy. If found, no update needed.
            return row
        except KeyError:
            pass
        new_tax_id = taxonomy._get_merged(current_tax_id)

        if new_tax_id and new_tax_id != current_tax_id:
            row['tax_id'] = new_tax_id
            log.warn('Replacing %s with %s [%s]', current_tax_id, new_tax_id,
                    row['seqname'])
        elif action == 'halt':
            raise KeyError("Unknown taxon {0}".format(current_tax_id))
        elif action == 'remove':
            logging.warn('Unknown Taxon: %s. Removing.', current_tax_id)
            row['tax_id'] = ''
        else:
            assert False

        return row

    return update_taxid


def action(args):
    """
    Show information about reference packages.
    """
    if not os.path.exists(args.database_file):
        log.error("Database does not exist: %s", args.database_file)
        return 1

    e = sqlalchemy.create_engine('sqlite:///{0}'.format(args.database_file))
    tax = Taxonomy(e, ncbi.ranks)

    headers, dialect, rows = load_csv(args.infile)

    for header in 'seqname', 'tax_id':
        if header not in headers:
            raise ValueError("Missing required field: {0}".format(header))

    update = taxid_updater(tax, args.unknown_action)
    updated = (update(row) for row in rows)

    with args.out_file as fp:
        writer = csv.DictWriter(fp, headers, dialect)
        writer.writeheader()
        writer.writerows(updated)

