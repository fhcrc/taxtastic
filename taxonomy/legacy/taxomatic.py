#!/usr/bin/env python

"""\
taxomatic.py
============

Usage: %prog action <options>

Creation, validation, and modification of reference packages for use
with `pplacer` and related software.

Actions
=======

Use `taxomatic.py -h` or `taxomatic.py --help` to print this help text.

help
  Print detailed help for an action below using `taxomatic.py help <action>`
create
  Create a reference package
check
  Verify that reference package is intact and valid.

"""

help_create = """\
action create
-------------

"""

help_check = """\
action check
------------

"""

help_tail = """\
Reference packages
==================

The full reference package specification can be seen at
http://github.com/fhcrc/taxtastic/wiki/refpkg.

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

import taxonomy
#from taxonomy.package import manifest_name, package_contents, write_config
from taxonomy import __version__ as version

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

    ## TODO: implement this block preceding `parser = OptionParser...`
    ## as a function or class elsewhere?

    valid_actions = ('help','create','check')

    short_usage = textwrap.dedent(__doc__)
    long_usage = short_usage + help_tail
    full_usage = short_usage + \
        '\n'.join(globals().get('help_'+x,'') for x in valid_actions) + help_tail

    argv = sys.argv[1:]

    args = []
    action = 'help'
    usage = short_usage
    print_opts = True

    if len(argv) == 0:
        pass
    elif argv[0] == 'help':
        usage = full_usage
        if len(argv) > 1:
            action, args = argv[1], []
            usage = short_usage + '\n' + globals().get('help_'+action,'')
    elif '-h' in argv or '--help' in argv:
        usage = long_usage
    else:
        action, args = argv[0], argv[1:]
        print_opts = False

    if action not in valid_actions:
        action_str = ','.join('"%s"'%a for a in valid_actions)
        msg = textwrap.fill(
            'Error: "%s" is not a valid action. Please choose from %s.' % \
                (action, action_str))
        sys.stderr.write(msg + '\n')
        sys.exit(1)

    parser = OptionParser(usage='',
                          version=version,
                          formatter=SimpleHelpFormatter())

    parser.set_defaults(
        package_name = './taxonomy.refpkg',
        verbose=0
        )

    parser.add_option("-v", "--verbose",
        action="count", dest="verbose",
        help="increase verbosity of screen output (eg, -v is verbose, -vv more so)")

    if action in ('create','check'):
        parser.add_option("-P", "--package-name",
            action="store", dest="package_name", type="string",
            help=xws("""Name of output directory. [default %default]
            """), metavar='PATH')

    if action in ('create'):
        parser.add_option("-t", "--tree-file",
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

        parser.add_option("-m", "--mask",
            action="store", dest="mask", type="string",
            help=xws("""
            Text file containing a mask.
            """), metavar='FILE')

    (options, args) = parser.parse_args(args=args)

    loglevel = {
        0:logging.WARNING,
        1:logging.INFO,
        2:logging.DEBUG
        }.get(options.verbose, logging.DEBUG)

    verbose_format = '%(levelname)s %(module)s %(lineno)s %(message)s'

    logformat = {0:'%(message)s',
        1:verbose_format,
        2:verbose_format}.get(options.verbose, verbose_format)

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    if action == 'help' or print_opts:
        print usage
        parser.print_help()
    else:
        try:
            if hasattr(taxonomy.package, action):
                getattr(taxonomy.package, action)(options)
            else:
                log.error('Sorry: the %s action is not yet implemented' % action)
                sys.exit(1)
        except OSError:
            log.error('A package named "%s" already exists' % options.package_name)
            sys.exit(2)

if __name__ == '__main__':
    main()
