#!/usr/bin/env python

"""\
============
taxomatic.py
============

Usage: taxomatic.py command <options>

This script accepts a set of files as input, performs some basic
sanity checks of format and content, and creates a package directory
with manifest.

Command line options
====================

Help text can be viewed using '-h' or '--help'.

"""

from optparse import OptionParser, IndentedHelpFormatter
import ConfigParser
import gettext
import logging
import os
import pprint
import sys
import textwrap
import time
import shutil
import hashlib

log = logging

import Taxonomy
from Taxonomy.package import manifest_name, package_contents, write_config

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
        package_name = './taxonomy.refpkg',
        verbose=0
        )

    parser.add_option("-P", "--package-name",
        action="store", dest="package_name", type="string",
        help=xws("""Name of output directory. [default %default]
        """), metavar='PATH')

    parser.add_option("--t", "--tree-file",
        action="store", dest="tree_file", type="string",
        help=xws("""
        Phylogenetic tree in newick format.
        """), metavar='FILE')

    parser.add_option("-s", "--tree-stats",
        action="store", dest="tree_stats", type="string",
        help=xws("""
        File containing tree statistics (for example "RAxML_info.whatever").
        """), metavar='FILE')

    parser.add_option("-f", "--aln-fasta",
        action="store", dest="aln_fasta", type="string",
        help=xws("""
        Multiple alignment in fasta format.
        """), metavar='FILE')

    parser.add_option("-S", "--aln-sto",
        action="store", dest="aln_sto", type="string",
        help=xws("""
        Multiple alignment in Stockholm format.
        """), metavar='FILE')

    parser.add_option("-p", "--profile",
        action="store", dest="profile", type="string",
        help=xws("""
        Alignment profile.
        """), metavar='FILE')

    parser.add_option("-H", "--hmm-profile",
        action="store", dest="hmm_profile", type="string",
        help=xws("""
        Alignment profile used by hmmer.
        """), metavar='FILE')

    parser.add_option("-i", "--seq-info",
        action="store", dest="seq_info", type="string",
        help=xws("""
        CSV format file decsribing the aligned reference sequences,
        minimally containing the fields "seqname" and "tax_id"
        """), metavar='FILE')

    parser.add_option("-T", "--taxonomy",
        action="store", dest="taxonomy", type="string",
        help=xws("""
        CSV format file defining the taxonomy. Fields include
        "tax_id","parent_id","rank","tax_name" followed by a column
        defining tax_id at each rank starting with root.
        """), metavar='FILE')

    parser.add_option("-v", "--verbose",
        action="count", dest="verbose",
        help="increase verbosity of screen output (eg, -v is verbose, -vv more so)")

    if not sys.argv[2:] or '-h' in sys.argv or '--help' in sys.argv:
        args = ['-h']
    else:
        command = sys.argv[1]
        args = sys.argv[2:]

    (options, args) = parser.parse_args(args=args)

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

    if command == 'create':
        try:
            Taxonomy.package.create(options.package_name, options)
        except OSError:
            log.error('A package named "%s" already exists' % options.package_name)
            sys.exit(2)

if __name__ == '__main__':
    main()
