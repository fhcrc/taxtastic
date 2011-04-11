#!/usr/bin/env python

"""\
taxtable.py
===========

Usage: %prog <options>

Create a csv file describing lineages for a set of taxa.

Command line options
--------------------

Help text can be viewed using the '-h' option.

"""

from optparse import OptionParser, IndentedHelpFormatter
import contextlib
import gettext
import logging
import os
import pprint
import sys
import textwrap

log = logging

try:
    from sqlalchemy import create_engine
except ImportError:
    create_engine = None
    print("""\n\n** Warning: this script requires the sqlalchemy package for some features; see "Installation." **\n""")
else:
    from sqlalchemy.exc import IntegrityError

import Taxonomy

class SimpleHelpFormatter(IndentedHelpFormatter):
    """Format help with indented section bodies.
    Modifies IndentedHelpFormatter to suppress leading "Usage:..."
    """
    def format_usage(self, usage):
        return gettext.gettext(usage)

def xws(s):
    return ' '.join(s.split())

def getlines(fname):
    with open(fname) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                yield line.strip()

def main():

    usage = textwrap.dedent(__doc__)

    parser = OptionParser(usage=usage,
                          version="$Id$",
                          formatter=SimpleHelpFormatter())

    parser.set_defaults(
        dest_dir = '.',
        dbfile = 'ncbi_taxonomy.db',
        new_database = False,
        source_name = 'unknown',
        verbose=0
        )

    parser.add_option("-o", "--outfile",
        action="store", dest="outfile", type="string",
        help=xws("""Output file containing lineages for the specified taxa (csv fomat);
        writes to stdout if unspecified"""), metavar='FILENAME')

    parser.add_option("-D", "--dest-dir",
        action="store", dest="dest_dir", type="string",
        help=xws("""Name of output directory.
        If --outfile is an absolute path, the path provided takes precedence for that
        file. [default is the current directory]"""), metavar='PATH')

    parser.add_option("-d", "--database-file",
        action="store", dest="dbfile", type="string",
        help=xws("""Filename of sqlite database [%default]."""), metavar='FILENAME')

    parser.add_option("-N", "--new-database", action='store_true',
                      dest="new_database", help=xws("""Include this
        option to overwrite an existing database and reload with
        taxonomic data (from the downloaded archive plus additional
        files if provided). [default %default]
        """))

    parser.add_option("-a", "--add-new-nodes", dest="new_nodes", help=xws("""
        An optional Excel (.xls) spreadsheet (requires xlrd) or
        csv-format file defining nodes to add to the
        taxonomy. Mandatory fields include
        "tax_id","parent_id","rank","tax_name"; optional fields
        include "source_name", "source_id". Other columns are ignored.
    """))

    parser.add_option("-S", "--source-name", dest="source_name", help=xws("""
        Names the source for new nodes. [default %default]
    """))

    parser.add_option("-t", "--tax-ids", dest="taxids", help=xws("""
        Specifies the taxids to include. Taxids may be provided in one
        of two ways: 1) as a comma delimited list on the command line
        (eg, '-t taxid1,taxid2'); or 2) as the name of a file
        containing a whitespace-delimited list of tax_ids (ie,
        separated by tabs, spaces, or newlines; lines beginning with
        "#" are ignored). May be omitted if --tax-names is used
        instead.
    """))

    parser.add_option("-n", "--tax-names", dest="taxnames", help=xws("""
        An optional file identifing taxa in the form of taxonomic
        names. Names are matched against both primary names and
        synonyms. Lines beginning with "#" are ignored. Taxa
        identified here will be added to those specified using
        --tax-ids.
    """))

    parser.add_option("-v", "--verbose",
        action="count", dest="verbose",
        help="increase verbosity of screen output (eg, -v is verbose, -vv more so)")

    parser.add_option("-q", "--quiet", action="store_true", dest="quiet",
        help="suppress output destined for stderr")


    (options, args) = parser.parse_args()

    loglevel = {
        0:logging.WARNING,
        1:logging.INFO,
        2:logging.DEBUG
        }.get(options.verbose, logging.DEBUG)

    verbose_format = '# %(levelname)s %(module)s %(lineno)s %(message)s'

    logformat = {0:'# %(message)s',
        1:verbose_format,
        2:verbose_format}.get(options.verbose, verbose_format)

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    # suppress output to stderr if --quiet specified
    if options.quiet:
        # Adapted from Alex Martelli's example on stackoverflow.com
        @contextlib.contextmanager
        def suppress_stderr():
            saved = sys.stderr
            class DevNull(object):
                def write(self, _): pass
            sys.stderr = DevNull()
            yield
            sys.stderr = saved
    
        sys.stderr = suppress_stderr()

    pth, fname = os.path.split(options.dbfile)
    dbname = options.dbfile if pth else os.path.join(options.dest_dir, fname)

    if not os.access(dbname, os.F_OK) or options.new_database:
        zfile = Taxonomy.ncbi.fetch_data(dest_dir=options.dest_dir, new=False)
        log.warning('creating new database in %s using data in %s' % (dbname, zfile))
        con = Taxonomy.ncbi.db_connect(dbname, new=True)
        Taxonomy.ncbi.db_load(con, zfile)
        con.close()
    else:
        log.warning('using taxonomy defined in %s' % dbname)

    if not create_engine:
        sys.exit('sqlalchemy is required, exiting.')

    engine = create_engine('sqlite:///%s' % dbname, echo = options.verbose > 1)
    tax = Taxonomy.Taxonomy(engine, Taxonomy.ncbi.ranks)

    # add nodes if necessary
    if options.new_nodes:
        log.warning('adding new nodes')
        new_nodes = Taxonomy.utils.get_new_nodes(options.new_nodes)
        for d in new_nodes:
            if options.source_name:
                d['source_name'] = options.source_name
                try:
                    tax.add_node(**d)
                except IntegrityError:
                    log.info('node with tax_id %(tax_id)s already exists' % d)

    # get a list of taxa
    taxa = set()

    taxids = options.taxids
    if taxids:
        if os.access(taxids, os.F_OK):
            log.warning('reading tax_ids from %s' % taxids)
            for line in getlines(taxids):
                taxa.update(set(line.split()))
        else:
            taxa = set([x.strip() for x in taxids.split(',')])

    taxnames = options.taxnames
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

    if options.outfile:
        pth, fname = os.path.split(options.outfile)
        csvname = options.outfile if pth else os.path.join(options.dest_dir, fname)
        log.warning('writing %s' % csvname)
        csvfile = open(csvname, 'w')
    else:
        csvfile = sys.stdout

    tax.write_table(None, csvfile = csvfile)

    engine.dispose()

if __name__ == '__main__':
    main()
