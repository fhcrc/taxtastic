#!/usr/bin/env python
import sys, os, string, argparse, re

# Insert one level above project directory to path for testing.
sys.path.insert(0, "../..")
from taxonomy.alignment import Alignment
#from Taxonomy.alignment import Alignment


def main():
    """
    Entry point for this script.
    """
    # Get command-line arguments.
    arguments = parse_arguments()
    min_length = arguments.min_length
    out_prefix = arguments.out_prefix
    reference_package = arguments.refpkg[0]
    sequence_files = arguments.seqfiles
    profile_version = arguments.profile_version
    # --squeeze is implicit when --mask is specified.
    #squeeze = arguments.squeeze or arguments.mask_file
    #mask = arguments.mask
    
    # Create alignment with hmmer for all sequence files.  Squeeze and 
    # use mask if desired.
    align = Alignment(reference_package=reference_package, 
                      out_prefix=out_prefix,
                      profile_version=profile_version,
                      min_length=min_length,
                     )

    align.hmmer_align(sequence_files=sequence_files, 
                      frag=True, ref=True,
                     )
                      #squeeze=squeeze, 


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='refpkg_align - align sequence ' + \
                                                 'files to be incorporated in a reference package')
    parser.add_argument('-o', '--outprefix', dest='out_prefix', help='Output file prefix. ' + \
                        'Defaults to refpkg_prefix.sequence_file_prefix.  Currently only works ' + \
                        'with a single sequence file')
    parser.add_argument('--profileversion', dest='profile_version', default=r"\d+", 
                        help='Integer or regular expression matching a profile version or version. ' + \
                        'For hmmer 3, "3" would be sufficient input.')
    parser.add_argument('--min-length', dest='min_length', type=int, default=1, metavar='N',
                        help='minimum sequence length. Defaults to 1.')
    parser.add_argument('refpkg', nargs=1, type=reference_package, help='Reference package directory')
    parser.add_argument('seqfiles', nargs='+', help='A list of one or more fasta files')
    # squeezing always done now
    #parser.add_argument('--squeeze', action='store_true', default=False, help='Squeeze out .s from an alignment ')
    # Mask done automatically if found in CONTENTS.json
    #parser.add_argument('--mask', dest='mask', action='store_true',  default=False,
    #                    help='Use mask specified in CONTENTS.json.  Implies --squeeze')

    arguments = parser.parse_args()
    # out_prefix only works when a single sequence file is passed. 
    if arguments.out_prefix and len(arguments.seqfiles) > 1:
        raise Exception, "--out-prefix cannot be used when more than a single sequence file is specified"

    return arguments



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
