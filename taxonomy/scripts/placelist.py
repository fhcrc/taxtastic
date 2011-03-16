#!/usr/bin/env python

"""\
placelist.py
============

Usage: %prog placefile <options>

Command line options
--------------------

Help text can be viewed using the '-h' option.

Note that the .place file is typically provided as the first argument,
but the '-f' option is provided for compatibility with frameworks
prohibiting positional arguments.
"""

from optparse import OptionParser, IndentedHelpFormatter
import csv
import gettext
import itertools
import json
import logging
import os
import pprint
import sys
import textwrap

log = logging

class SimpleHelpFormatter(IndentedHelpFormatter):
    """Format help with indented section bodies.
    Modifies IndentedHelpFormatter to suppress leading "Usage:..."
    """
    def format_usage(self, usage):
        return gettext.gettext(usage)

def xws(s):
    return ' '.join(s.split())

# def read(fname, ncol):

#     counter = itertools.count

#     with open(fname) as lines:
#         data = []
#         skip = False
#         for line in lines:
#             if line.startswith('#'):
#                 continue
#             elif skip:
#                 skip = False
#             elif line.startswith('>'):
#                 if data:
#                     yield data
#                 seqname = line.strip('\n >')
#                 skip = True ## skip the next line (the sequence string)
#                 data = []
#                 count = counter(1)
#             else:
#                 spl = line.split()[:-1]                
#                 tax_id = spl[-1]
#                 spl[-1] = tax_id.split('*')[1] if tax_id != 'none' else None
#                 data.append([seqname, str(count.next())] + spl)

#         yield data

def read(fname):

    """
    Returns an iterator over pplacer results.
    """

    with open(fname) as infile:
        data = json.load(infile)

    placements = data['placements']

    yield ['name','hit'] + data['fields']    
    for d in placements:
        name = d['n']
        hit = itertools.count(1)
        for p in d['p']:
            yield name + [hit.next()] + p

            
def main():

    usage = textwrap.dedent(__doc__)

    parser = OptionParser(usage=usage,
                          version="$Id$",
                          formatter=SimpleHelpFormatter())

    parser.set_defaults(
        infile = None,
        outfile = None
        )

    parser.add_option("-f", "--infile",
        action="store", dest="infile", type="string",
        help=xws("""Input .place file"""))

    parser.add_option("-o", "--outfile",
        action="store", dest="outfile", type="string",
        help=xws("""Output CSV file containing pplacer data."""))

    parser.add_option("-v", "--verbose",
        action="count", dest="verbose",
        help="increase verbosity of screen output (eg, -v is verbose, -vv more so)")

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

    try:
        fname = args[0]
    except IndexError:
        fname = options.infile

    if not fname:
        print('Error: an input file is required.\n')
        parser.print_usage()
        exit(1)
        
    fout = open(options.outfile, 'w') if options.outfile else sys.stdout
    writer = csv.writer(fout, quoting=csv.QUOTE_NONNUMERIC, delimiter=',')
    writer.writerows(read(fname))
    
if __name__ == '__main__':
    main()
