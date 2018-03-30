#!/usr/bin/env python

"""Generate help text for each subcommand.

Output files are saved in ``outdir/{subcommand}.txt``
"""

import os
import sys
import argparse
import subprocess

from taxtastic import subcommands


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('taxit', help='path to taxit executable')
    parser.add_argument('-d', '--outdir', help="Output directory")

    args = parser.parse_args(arguments)

    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

    # execute each subcommand in a subshell because calling
    # taxit.main(args) exits after printing help text for first
    # subcommand.
    for name, module in subcommands.itermodules():
        with open(os.path.join(args.outdir, name + '.txt'), 'w') as f:
            subprocess.call([args.taxit, name, '-h'], stdout=f)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
