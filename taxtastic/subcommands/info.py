"""
Show information about reference packages.
"""
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

import logging
import os.path
import shutil

from Bio import Phylo, SeqIO

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-n', '--seq-names', action = 'store_true',
                        help = 'print a list of sequence names')

def tree_names(pkg):
    tree = Phylo.read(pkg.file_abspath('tree'), 'newick')    
    return [x.name for x in tree.get_terminals()]
    
def action(args):
    """
    Show information about reference packages.
    """
    log.info('loading reference package')

    pkg = refpkg.Refpkg(args.refpkg)    
    tnames = tree_names(pkg)

    with open(pkg.file_abspath('aln_sto')) as f:
        seqs = SeqIO.parse(f, 'stockholm')
        snames = [s.id for s in seqs]

    if not set(tnames) == set(snames):
        sys.exit('Error: sequence names in the Stockholm alignment and tree file differ.')
        
    if args.seq_names:    
        for name in tnames:
            print name
