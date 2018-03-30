#!/usr/bin/env python
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
Creation, validation, and modification of reference packages for use
with `pplacer` and related software.
"""

import argparse
from argparse import RawDescriptionHelpFormatter
import sys
import os
import logging
from taxtastic import subcommands, __version__ as version

DESCRIPTION = __doc__.strip()


def main(argv=None):
    argv = argv or sys.argv[1:]
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
    logging.basicConfig(stream=sys.stderr, format=logformat, level=loglevel)

    try:
        return action(arguments)
    finally:
        subcommands.close_all_files(arguments)


def parse_arguments(argv):
    """Create the argument parser

    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    base_parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('-V', '--version', action='version',
                        version='taxit v' + version,
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

    for name, mod in subcommands.itermodules(
            os.path.split(subcommands.__file__)[0]):
        # set up subcommand help text. The first line of the dosctring
        # in the module is displayed as the help text in the
        # script-level help message (`script -h`). The entire
        # docstring is displayed in the help message for the
        # individual subcommand ((`script action -h`)).
        subparser = subparsers.add_parser(
            name,
            prog='taxit {}'.format(name),
            help=mod.__doc__.lstrip().split('\n', 1)[0],
            description=mod.__doc__,
            formatter_class=RawDescriptionHelpFormatter,
            parents=[base_parser])

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
