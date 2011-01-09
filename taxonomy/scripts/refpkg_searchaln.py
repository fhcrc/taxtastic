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
    sequence_file = arguments.seqfile[0]
    search_options = arguments.search_options    

    align = Alignment(reference_package=reference_package, 
                      out_prefix=out_prefix,
                     )


    hmmsearch_output_file = [ align.hmmer_search(search_options=search_options, sequence_file=sequence_file) ]
 
    # Create alignment with hmmer for the file containing 
    # recruited sequences.  Squeeze and use mask if desired.
    align.hmmer_align(sequence_files=hmmsearch_output_file,
                      squeeze=True, mask=True, 
                      frag=True, ref=False,
                      separate_steps=True,
                      sequence_file_format='stockholm',
                     )


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='refpkg_searchaln.py - hmmsearch wrapper ' + \
                                                 'with refpkg-align.py-like functionality.')
    parser.add_argument('-o', '--outprefix', dest='out_prefix', help='Output file prefix. ' + \
                        'Defaults to refpkg_prefix.sequence_file_prefix.  Currently only works ' + \
                        'with a single sequence file')
    parser.add_argument('--search-opts', dest='search_options', metavar='OPTS', help='hmmsearch options, such as "-E 1"')
    parser.add_argument('refpkg', nargs=1, type=reference_package, help='Reference package directory')
    parser.add_argument('seqfile', nargs=1, help='A single fasta files')
    return parser.parse_args()


# This argparse type was copied from refpkg_align.py, so it would probably make more
# sense to put it somewhere else...
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
