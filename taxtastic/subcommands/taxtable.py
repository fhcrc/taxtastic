"""Creates a CSV file describing lineages for a set of taxa"""

from taxtastic import ncbi
from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import getlines, get_new_nodes

from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
import os.path
import sys
from textwrap import dedent

import logging
log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('-a', '--add-new-nodes', dest='new_nodes',
        help='An optional Excel (.xls) spreadsheet (requires xlrd) or '
             'csv-format file defining nodes to add to the taxonomy. '
             'Mandatory fields include "tax_id","parent_id","rank",'
             '"tax_name"; optional fields include "source_name", '
             '"source_id". Other columns are ignored.')

    parser.add_argument('-d', '--database-file',
        action='store', dest='dbfile', default='ncbi_taxonomy.db',
        help='Filename of sqlite database [%(default)s].',
        metavar='FILENAME')

    parser.add_argument('-D', '--dest-dir', action='store',
        default='.', dest='dest_dir', metavar='PATH',
        help='Name of output directory. If --out-file '
             'is an absolute path, the path provided takes precedence '
             'for that file. [default is the current directory]')

    parser.add_argument('-o', '--out-file',
        action='store', dest='out_file', metavar='FILENAME',
        help='Output file containing lineages for the specified taxa '
             '(csv fomat); writes to stdout if unspecified')

    parser.add_argument('-n', '--tax-names', dest='taxnames',
        help='An optional file identifing taxa in the form of taxonomic '
             'names. Names are matched against both primary names and '
             'synonyms. Lines beginning with "#" are ignored. Taxa '
             'identified here will be added to those specified using '
             '--tax-ids')

    parser.add_argument('-N', '--new-database', action='store_true',
        dest='new_database', default=False, help='Include this '
             'option to overwrite an existing database and reload with '
             'taxonomic data (from the downloaded archive plus additional '
             'files if provided). [default %(default)s]')
        
    parser.add_argument('-S', '--source-name', dest='source_name',
        default='unknown', help='Names the source for new nodes. '
                '[default %(default)s]')

    parser.add_argument('-t', '--tax-ids', dest='taxids',
        help='Specifies the taxids to include. Taxids may be '
             'provided in one of two ways: 1) as a comma delimited list '
             'on the command line (eg, "-t taxid1,taxid2"); or 2) as the '
             'name of a file containing a whitespace-delimited list of '
             'tax_ids (ie, separated by tabs, spaces, or newlines; lines '
             'beginning with "#" are ignored). May be omitted if '
             '--tax-names is used instead')

    parser.add_argument('-z', '--clobber-zipfile', action='store_true',
                        dest='clobber_zipfile', default=False,
                        help = dedent("""Download a new zip archive
                        containing NCBI taxonomy even if one already
                        exists [default %(default)s]""")
                        )
    
def action(arguments):
    pth, fname = os.path.split(arguments.dbfile)
    dbname = arguments.dbfile if pth else os.path.join(arguments.dest_dir, fname)

    zfile, downloaded = ncbi.fetch_data(
        dest_dir = arguments.dest_dir,
        clobber = arguments.clobber_zipfile)

    if not os.access(dbname, os.F_OK) or arguments.new_database:
        log.warning('creating new database in %s using data in %s' % (dbname, zfile))
        con = ncbi.db_connect(dbname, clobber=True)
        ncbi.db_load(con, zfile)
        con.close()
    else:
        log.warning('using taxonomy defined in %s' % dbname)

    sys.exit()
        
    engine = create_engine('sqlite:///%s' % dbname, echo=arguments.verbosity > 2)
    tax = Taxonomy(engine, ncbi.ranks)

    # add nodes if necessary
    if arguments.new_nodes:
        log.warning('adding new nodes')
        new_nodes = get_new_nodes(arguments.new_nodes)
        for d in new_nodes:
            if arguments.source_name:
                d['source_name'] = arguments.source_name
                try:
                    tax.add_node(**d)
                except IntegrityError:
                    log.info('node with tax_id %(tax_id)s already exists' % d)

    # get a list of taxa
    taxa = set()

    taxids = arguments.taxids
    if taxids:
        if os.access(taxids, os.F_OK):
            log.warning('reading tax_ids from %s' % taxids)
            for line in getlines(taxids):
                taxa.update(set(line.split()))
        else:
            taxa = set([x.strip() for x in taxids.split(',')])

    taxnames = arguments.taxnames
    if taxnames:
        for tax_name in getlines(taxnames):
            tax_id, primary_name, is_primary = tax.primary_from_name(tax_name)
            taxa.add(tax_id)
            if not is_primary:
                log.warning('%(tax_id)8s  %(tax_name)40s -(primary name)-> %(primary_name)s' % locals())

    log.warning('calculating lineages for %s taxa' % len(taxa))
    for taxid in taxa:
        log.warning('adding %s' % taxid)

        tax.lineage(taxid)

    if arguments.out_file:
        pth, fname = os.path.split(arguments.out_file)
        csvname = arguments.out_file if pth else os.path.join(arguments.dest_dir, fname)
        log.warning('writing %s' % csvname)
        csvfile = open(csvname, 'w')
    else:
        csvfile = sys.stdout

    tax.write_table(None, csvfile = csvfile)

    engine.dispose()
