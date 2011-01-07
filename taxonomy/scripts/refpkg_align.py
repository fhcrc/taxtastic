#!/usr/bin/env python
import sys, os, string, argparse, re

# Insert one level above project directory to path for testing.
#sys.path.insert(0, "../..")
#from taxonomy.alignment import Alignment
from Taxonomy.alignment import Alignment


def main():
    """
    Entry point for this script.
    """
    # Get command-line arguments.
    arguments = parse_arguments()
    out_prefix = arguments.out_prefix
    reference_package = arguments.refpkg[0]
    sequence_files = arguments.seqfiles
    # --squeeze is implicit when --mask is specified.
    squeeze = arguments.squeeze or arguments.mask_file
    mask = arguments.mask
    
    # Create alignment with hmmer for all sequence files.  Squeeze and 
    # use mask if desired.
    align = Alignment(reference_package=reference_package, 
                      out_prefix=out_prefix,
                     )

    align.hmmer_align(sequence_files=sequence_files, 
                      squeeze=squeeze, mask=mask, 
                      frag=True, ref=True,
                     )


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='refpkg_align - align sequence ' + \
                                                 'files to be incorporated in a reference package')
    parser.add_argument('-o', '--outprefix', dest='out_prefix', help='Output file prefix. ' + \
                        'Defaults to refpkg_prefix.sequence_file_prefix.  Currently only works ' + \
                        'with a single sequence file')
    parser.add_argument('refpkg', nargs=1, type=reference_package, help='Reference package directory')
    parser.add_argument('seqfiles', nargs='+', help='A list of one or more fasta files')
    parser.add_argument('--squeeze', action='store_true', default=False, help='Squeeze out .s from an alignment ')
    parser.add_argument('--mask', dest='mask', action='store_true',  default=False,
                        help='Use mask specified in CONTENTS.json.  Implies --squeeze')
    return parser.parse_args()


def reference_package(reference_package):
    """
    A custom argparse 'type' to make sure the path to the reference package exists.
    """
    if os.path.isdir(reference_package):
        # Remove trailing / from directory if found.
        if reference_package.endswith('/'):
            reference_package = reference_package[:-1]
        return reference_package
    else:
        raise Exception, 'Path to reference package does not exist: ' + reference_package

if __name__ == '__main__':
    sys.exit(main())
