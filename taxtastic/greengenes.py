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
"""
Methods and variables specific to the GreenGenes taxonomy.
"""
import contextlib
import logging
import itertools
import sqlite3
import tarfile

from . import taxdb
from .errors import IntegrityError

log = logging

ranks = taxdb.ranks + ['otu']

# Map from rank abbreviation to taxdb rank name
_rank_map = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
        'f': 'family', 'g': 'genus', 's': 'species'}

greengenes_data_url = 'http://greengenes.lbl.gov/Download/Sequence_Data/' \
        'Fasta_data_files/Caporaso_Reference_OTUs/gg_otus_4feb2011.tgz'
tax_map = 'gg_otus_4feb2011/taxonomies/greengenes_tax.txt'

db_connect = taxdb.db_connect

def _get_source_id(con):
    """Get the source.id corresponding to greengenes"""
    cursor = con.cursor()
    cursor.execute('SELECT id FROM source where name = ?', ['greengenes'])
    return cursor.fetchone()

def _load_rows(rows, con):
    """
    Load GreenGenes taxonomy into database.

    Sequences are labeled with tax_ids starting from GG00000001

    No action is taken if tables are already populated.
    """
    source_id = _get_source_id(con)
    cursor = con.cursor()
    if taxdb.has_row(cursor, 'nodes') or taxdb.has_row(cursor, 'names'):
        log.warning("nodes table has data. not updating.")
        return

    count = itertools.count()

    def _get_or_insert(name, rank, parent_id, source_id):
        # Check for existence
        cursor.execute("SELECT tax_id from names WHERE tax_name = ?", [name])
        tax_id = cursor.fetchone()
        if tax_id:
            return tax_id[0]

        # No tax_id found - create
        tax_id = 'GG{0:08d}'.format(count.next())

        cursor.execute("""INSERT INTO nodes (tax_id, parent_id, rank, source_id)
                VALUES (?, ?, ?, 1)""", (tax_id, parent_id, rank))
        cursor.execute("""INSERT INTO names (tax_id, tax_name, is_primary)
                VALUES (?, ?, 1)""", (tax_id, name))
        return tax_id

    with con:
        for classes in rows:
            parent = None
            for rank, name in classes:
                parent = _get_or_insert(name, rank, parent, source_id)

def db_load(con, archive, maxrows=None):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.
    """
    try:
        with tar_member(archive) as handle:
            rows = _parse_gg(handle)
            _load_rows(rows, con)
    except sqlite3.IntegrityError, err:
        raise IntegrityError(err)

def fetch_data(dest_dir='.', clobber=False, url=greengenes_data_url):
    """
    Download data from greengenes required to generate local taxonomy
    database.

    * dest_dir - directory in which to save output files (created if necessary).
    * clobber - don't download if False and target of url exists in dest_dir
    * url - url to archive;

    Returns (fname, downloaded), where fname is the name of the
    downloaded zip archive, and downloaded is True if a new files was
    downloaded, false otherwise.
    """
    return taxdb.fetch_url(url, dest_dir, clobber)

def _parse_classes(classes, otu_id):
    """
    Parse classes from GreenGenes taxonomy

    classes - split list of GreenGenes taxonomy assignments, e.g.
       ['k__Bacteria', 'p__Proteobacteria']
    otu_id - OTU ID

    Returns list of (rank, name) tuples
    """
    split = [i.split('__') for i in classes]
    result = [(_rank_map[cls_key], cls or None)
             for i, (cls_key, cls) in enumerate(split) if cls]

    # special handling for species -> genus
    if result[-1][0] == 'species' and result[-2][0] == 'genus':
        species = result[-1][1]
        genus = result[-2][1]
        if species.startswith(genus):
            l = len(genus)
            species = ' '.join((species[:l], species[l:]))

    result.append(('otu', str(otu_id)))
    return result

def _parse_gg(handle):
    """
    Parse GreenGenes data

    Yields parsed classes
    """
    for line in handle:
        otu_id, classes = line.rstrip().split('\t')
        otu_id = int(otu_id)
        classes = _parse_classes(classes.split(';'), otu_id)
        yield classes

@contextlib.contextmanager
def tar_member(archive, fname=tax_map):
    """
    Provides a file-like object for the path fname in tar archive
    """
    with tarfile.open(archive) as tfile:
        yield tfile.extractfile(fname)
