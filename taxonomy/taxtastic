#!/usr/bin/env python

"""\
taxtastic
=========

Usage: %prog action <options>

Creation, validation, and modification of reference packages for use
with `pplacer` and related software.

Use `taxomatic.py -h` or `taxomatic.py --help` to print help text.
"""

import argparse
import sys
import os
import logging
from taxonomy import package, __version__ as version

log = logging
PROG = os.path.basename(__file__)
DESCRIPTION = 'To Be Named -- Creation, validation, and modification of ' + \
              'reference packages for use with `pplacer` and related software.'

# Noah: do we need to do this check or is the sqlalchemy requirement in 
# setup.py sufficient?
try:
    from sqlalchemy import create_engine
except ImportError:
    create_engine = None
    print("\n\nWarning: this script requires the sqlalchemy package for " + \
          "some features; see \"Installation section in the taxtastic " + \
          "documentation for more details.\" \n")
else:
    from sqlalchemy.exc import IntegrityError


def main():
    """
    """
    action, arguments = parse_arguments()

    loglevel = {
        0:logging.WARNING,
        1:logging.INFO,
        2:logging.DEBUG
        }.get(arguments.verbose, logging.DEBUG)

    verbose_format = '%(levelname)s %(module)s %(lineno)s %(message)s'

    logformat = {0:'%(message)s',
        1:verbose_format,
        2:verbose_format}.get(arguments.verbose, verbose_format)

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    try:
        if hasattr(package, action):
            getattr(package, action)(arguments)
        else:
            log.error('Sorry: the %s action is not yet implemented' % action)
            sys.exit(1)
    except OSError:
        log.error('A package named "%s" already exists' % arguments.package_name)
        sys.exit(2)


def parse_arguments(action_arguments=None):
    """
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description=DESCRIPTION, prog=PROG)

    parser.add_argument('-V', '--version', action='version',
            version='%(prog)s v' + version,
            help='Print the version number and exit')

    ##########################
    # Setup all sub-commands #
    ########################## 

    subparsers = parser.add_subparsers(dest='subparser_name')
    
    # Begin help sub-command 
    parser_help = subparsers.add_parser('help', help='Detailed help for ' + \
        'actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)
    # End help sub-command

    # Begin create sub-command
    parser_create = subparsers.add_parser('create', 
        help='Create a reference package')

    parser_create.add_argument("-a", "--author",
        action="store", dest="author",
        help='Person who created the reference package', metavar='NAME')

    parser_create.add_argument("-d", "--description",
        action="store", dest="description",
        help='An arbitrary description field', metavar='TEXT')

    parser_create.add_argument("-f", "--aln-fasta",
        action="store", dest="aln_fasta",
        help='Multiple alignment in fasta format', metavar='FILE')

    parser_create.add_argument("-i", "--seq-info",
        action="store", dest="seq_info",
        help='CSV format file describing the aligned reference ' + \
             'sequences, minimally containing the fields "seqname" ' + \
             'and "tax_id"', metavar='FILE')

    parser_create.add_argument("-l", "--locus",
        action="store", dest="locus", required=True,
        help='The locus described by the reference package', metavar='LOCUS')

    parser_create.add_argument("-m", "--mask",
        action="store", dest="mask",
        help='Text file containing a mask', metavar='FILE')

    parser_create.add_argument("-p", "--profile",
        action="store", dest="profile",
        help='Alignment profile', metavar='FILE')

    parser_create.add_argument('-P', '--package-name', 
        action='store', dest='package_name',
        default='./taxonomy.refpkg', metavar='PATH',
        help='Name of output directory [default %(default)s]')

    parser_create.add_argument("-r", "--package-version",
        action="store", dest="package_version",
        help='Release version for the reference package', metavar='VERSION')

    parser_create.add_argument("-s", "--tree-stats",
        action="store", dest="tree_stats",
        help='File containing tree statistics (for example ' + \
             'RAxML_info.whatever")', metavar='FILE')

    parser_create.add_argument("-S", "--aln-sto",
        action="store", dest="aln_sto",
        help='Multiple alignment in Stockholm format', metavar='FILE')

    parser_create.add_argument("-t", "--tree-file",
        action="store", dest="tree_file",
        help='Phylogenetic tree in newick format', 
        metavar='FILE')

    parser_create.add_argument("-T", "--taxonomy",
        action="store", dest="taxonomy",
        help='CSV format file defining the taxonomy. Fields include ' + \
             '"tax_id","parent_id","rank","tax_name" followed by a column ' + \
             'defining tax_id at each rank starting with root', metavar='FILE')
    # End create sub-command

    # Begin check sub-command
    parser_check = subparsers.add_parser('check', 
        help='The check action is not yet implemented')
        #help='Verify that a reference package is intact and valid')
    parser_check.add_argument('-P', '--package-name', 
        action='store', dest='package_name',
        default='./taxonomy.refpkg', metavar='PATH',
        help='Name of output directory. [default %(default)s]')
    # End check sub-command

    # Begin taxtable sub-command
    parser_help = subparsers.add_parser('taxtable', help='Create a CSV ' + \
        'file describing lineages for a set of taxa')

    parser.add_argument('-a', '--add-new-nodes', dest='new_nodes', 
        help='An optional Excel (.xls) spreadsheet (requires xlrd) or ' + \
             'csv-format file defining nodes to add to the taxonomy. ' + \
             'Mandatory fields include "tax_id","parent_id","rank",' + \
             '"tax_name"; optional fields include "source_name", ' + \
             '"source_id". Other columns are ignored.')

    parser.add_argument('-d', '--database-file',
        action='store', dest='dbfile', default='ncbi_taxonomy.db',
        help='Filename of sqlite database [%(default)s].'), 
        metavar='FILENAME')

    parser.add_argument('-D', '--dest-dir', action='store', 
        dest='dest_dir', help='Name of output directory. If --out-file ' + \
             'is an absolute path, the path provided takes precedence ' + \
             'for that file. [default is the current directory]', 
        metavar='PATH')

    parser.add_argument('-o', '--out-file',
        action='store', dest='out_file', type='string',
        help='Output file containing lineages for the specified taxa ' + \
             '(csv fomat); writes to stdout if unspecified', 
        metavar='FILENAME')

    parser.add_argument('-N', '--new-database', action='store_true',
        dest='new_database', default=False, help='Include this ' + \
             'option to overwrite an existing database and reload with ' + \
             'taxonomic data (from the downloaded archive plus additional ' + \
             'files if provided). [default %(default)s]')

    parser.add_argument('-S', '--source-name', dest='source_name', 
        default='unknown', help='Names the source for new nodes. ' + \
                '[default %(default)s]')

    parser.add_argument('-t', '--tax-ids', dest='taxids', 
        help='Specifies the taxids to include. Taxids may be ' + \
             'provided in one of two ways: 1) as a comma delimited list ' + \
             'on the command line (eg, '-t taxid1,taxid2'); or 2) as the ' + \
             'name of a file containing a whitespace-delimited list of ' + \
             'tax_ids (ie, separated by tabs, spaces, or newlines; lines ' + \
             'beginning with "#" are ignored). May be omitted if ' + \
             '--tax-names is used instead.')

    parser.add_argument("-n", "--tax-names", dest="taxnames", help=xws("""
        An optional file identifing taxa in the form of taxonomic
        names. Names are matched against both primary names and
        synonyms. Lines beginning with "#" are ignored. Taxa
        identified here will be added to those specified using
        --tax-ids.
    """))



    # End taxtable sub-command

    # With the exception of 'help', all subcommands can share a  
    # number of arguments, which are added here.
    for subcommand in subparsers.choices.keys():
        if subcommand == 'help': continue
        subparser = subparsers.choices[subcommand]
        subparser.add_argument('-v', '--verbose', action='count', dest='verbose', 
            help='Increase verbosity of screen output (eg, -v is verbose,' + \
                 '-vv more so)', default=0)
        parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
            help="Suppress output destined for stderr")


    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or* 
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(action_arguments) if parser.parse_args(action_arguments) else parser.parse_args()
    # Determine which action is in play.
    action = arguments.subparser_name

    # Support help <action> by simply having this function call itself and 
    # translate the arguments into something that argparse can work with.
    if action == 'help':
        parse_arguments(action_arguments=[str(arguments.action[0]), '-h'])

    return action, arguments



if __name__ == '__main__':
    sys.exit(main())
