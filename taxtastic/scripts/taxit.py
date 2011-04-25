#!/usr/bin/env python
"""\
taxtastic.py
============

Usage: %prog action <options>

Creation, validation, and modification of reference packages for use
with `pplacer` and related software.  Can also create CSV files
describing lineages for a set of taxa.

Use `taxomatic.py -h` or `taxomatic.py --help` to print help text.
"""

import argparse
import sys
import os
import logging
from taxtastic import subcommands, __version__ as version

PROG = os.path.basename(__file__)
DESCRIPTION = 'taxtastic.py -- Creation, validation, and modification of ' + \
              'reference packages for use with `pplacer` and related software.'

def main(argv):
    action, arguments = parse_arguments(argv)

    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(arguments.verbosity, logging.DEBUG)

    if arguments.verbosity > 1:
        logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
    else:
        logformat = '%(message)s'

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    return action(arguments)

def parse_arguments(argv):
    """
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description=DESCRIPTION, prog=PROG)

    parser.add_argument('-V', '--version', action='version',
        version='%(prog)s v' + version,
        help='Print the version number and exit')

    parser.add_argument('-v', '--verbose',
        action='count', dest='verbosity', default=1,
        help='Increase verbosity of screen output (eg, -v is verbose, '
             '-vv more so)')
    parser.add_argument('-q', '--quiet',
        action='store_const', dest='verbosity', const=0,
        help='Suppress output')


    ##########################
    # Setup all sub-commands #
    ##########################

    subparsers = parser.add_subparsers(dest='subparser_name')

    # Begin help sub-command
    parser_help = subparsers.add_parser(
        'help', help='Detailed help for actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)
    # End help sub-command

    actions = {}
    for name, mod in subcommands.itermodules():
        subparser = subparsers.add_parser(name, help=mod.__doc__)
        mod.build_parser(subparser)
        actions[name] = mod.action

    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or*
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(argv)
    # Determine which action is in play.
    action = arguments.subparser_name

    # Support help <action> by simply having this function call itself and
    # translate the arguments into something that argparse can work with.
    if action == 'help':
        return parse_arguments([str(arguments.action[0]), '-h'])

    return actions[action], arguments
