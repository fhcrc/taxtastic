"""Creates a reference package"""

import logging
import os
import time
import shutil
import hashlib
import re
import json
from collections import defaultdict

from taxtastic import utils
from taxtastic import package

log = logging.getLogger(__name__)

FORMAT_VERSION = '1.0'

MANIFEST_NAME = 'CONTENTS.json'
PHYLO_MODEL_FILE = 'phylo_model.json'

PACKAGE_CONTENTS = {
    'metadata': [
        'create_date', 'author', 'description', 'package_version',
        'empirical_frequencies', 'locus', 'format_version',
    ],
    'files': [
        'tree_file', 'tree_stats', 'aln_fasta', 'aln_sto',
        'profile', 'seq_info', 'taxonomy', 'mask', 'phylo_model_file',
    ],
    'md5': []
}

class ConfigError(Exception):
    pass

class StatsParser(object):
    """
    StatsParser class - parse tree stats output files generated
    by RAxML (2 types) and PhyML.
    """

    def __init__(self, file_name):
        """
        * file_name - name of input file to parse.

        Supported file types include:
        * raxml_condensed - subs_rates on a single line
        * raxml_re_estimated - subs_rates 1 per line
        * raxml_aa - amino acid
        * phyml_dna
        * phyml_aa
        """

        self.file_name = file_name
        self.file_type = 'unknown'  # Type of file examined, defaults to unknown.
        self.input_text = '' # Full text of input file.
        # Property of the class that gets populated when matches are found
        # within an output file.
        self.stats_values = defaultdict(dict)


### Public methods ###

    def get_stats_values(self):
        """
        Return contents of stats_values, in the format of a
        dictionary.
        """
        return self.stats_values


    def get_stats_json(self):
        """
        Return contents of stats_values in JSON format.
        """
        return json.dumps(self.stats_values, indent=2)


    def write_stats_json(self, out_file):
        """
        Write the values parsed from the statistics file in JSON
        format to out_file.
        """

        with open(out_file, 'w') as out:
            out.write(self.get_stats_json())


    # Return contents of file_type.
    def get_file_type(self):
        return self.file_type


    def parse_stats_data(self):
        """
        Extract data from stats output files.
        """

        # Read input file into a string for multiline matching.
        with open(self.file_name, 'r') as input_file:
            self.input_text = input_file.read()

        match_found = False
        for file_type in ['raxml_condensed','raxml_re_estimated','raxml_aa','phyml_dna','phyml_aa']:
            match_found = getattr(self, '_parse_'+file_type)()
            if match_found:
                self.file_type = file_type
                # Determine if empirical frequencies were used for the tree.
                self.stats_values['empirical_frequencies'] = self.is_empirical()
                break

        return match_found


    def is_empirical(self):
        """
        Determine if empirical base frequencies were used.  Defaults to true
        unless proven otherwise.  Only matches RaxML-generated files.
        """
        # Default to Ture for empirical_frequencies, unless proven otherwise.
        empirical_frequencies = True


        # Determine if the file type is RaxML.  If nucleotides, leave
        # empirical_frequencies true.  Otherwise look for the
        # string "Empirical Base Frequencies".  If it isn't there,
        # empirical_frequencies gets set to false.
        if self.file_type.startswith('raxml'):
            # empirical_frequencies will always be true for nucleotides.
            if self.stats_values['datatype'] != 'DNA' and \
               self.stats_values['datatype'] != 'RNA':

                regex = re.compile('.*(RAxML version .*?)Empirical Base Frequencies:',
                        re.M|re.DOTALL)

                if not regex.match(self.input_text):
                    empirical_frequencies = False
        else:
            print "Warning: phyml stats files don't specify empirical or " + \
                  "model frequencies; assuming empirical."


        return empirical_frequencies


### Private methods ###

    def _parse_raxml_condensed(self):
        """
        Parse RAxML info file - condensed.
        """

        regex = re.compile('.*(RAxML version .*?) release.*DataType: (\w+).*Substitution Matrix: (\w+).*alpha\[0\]: (\d+\.\d+) rates\[0\] ac ag at cg ct gt: (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) ',
                           re.M|re.DOTALL)

        if (regex.match(self.input_text)):
            self.stats_values['program'] = regex.match(self.input_text).group(1) # program
            self.stats_values['datatype'] = regex.match(self.input_text).group(2) # datatype
            self.stats_values['subs_model'] = regex.match(self.input_text).group(3) # subs_model
            self.stats_values['ras_model'] = 'gamma' # gamma is present
            self.stats_values['gamma']['alpha'] = float(regex.match(self.input_text).group(4)) # alpha
            self.stats_values['gamma']['n_cats'] = 4 # n_cats - default is 4 for RaxML
            self.stats_values['subs_rates']['ac'] = float(regex.match(self.input_text).group(5)) # ac
            self.stats_values['subs_rates']['ag'] = float(regex.match(self.input_text).group(6)) # ag
            self.stats_values['subs_rates']['at'] = float(regex.match(self.input_text).group(7)) # at
            self.stats_values['subs_rates']['cg'] = float(regex.match(self.input_text).group(8)) # cg
            self.stats_values['subs_rates']['ct'] = float(regex.match(self.input_text).group(9)) # ct
            self.stats_values['subs_rates']['gt'] = float(regex.match(self.input_text).group(10)) # gt
            return True
        else:
            return False

    def _parse_raxml_re_estimated(self):
        """
        Parse RAxML info file - re-estimated.
        """

        regex = re.compile('.*(RAxML version .*?) release.*DataType: (\w+).*Substitution Matrix: (\w+).*alpha: (\d+\.\d+).*rate A <-> C: (\d+\.\d+).*rate A <-> G: (\d+\.\d+).*rate A <-> T: (\d+\.\d+).*rate C <-> G: (\d+\.\d+).*rate C <-> T: (\d+\.\d+).*rate G <-> T: (\d+\.\d+)',
                           re.M|re.DOTALL)

        if (regex.match(self.input_text)):
            self.stats_values['program'] = regex.match(self.input_text).group(1) # program
            self.stats_values['datatype'] = regex.match(self.input_text).group(2) # datatype
            self.stats_values['subs_model'] = regex.match(self.input_text).group(3) # subs_model
            self.stats_values['ras_model'] = 'gamma' # gamma is present
            self.stats_values['gamma']['alpha'] = float(regex.match(self.input_text).group(4)) # alpha
            self.stats_values['gamma']['n_cats'] = 4 # n_cats - default is 4 for RaxML
            self.stats_values['subs_rates']['ac'] = float(regex.match(self.input_text).group(5)) # ac
            self.stats_values['subs_rates']['ag'] = float(regex.match(self.input_text).group(6)) # ag
            self.stats_values['subs_rates']['at'] = float(regex.match(self.input_text).group(7)) # at
            self.stats_values['subs_rates']['cg'] = float(regex.match(self.input_text).group(8)) # cg
            self.stats_values['subs_rates']['ct'] = float(regex.match(self.input_text).group(9)) # ct
            self.stats_values['subs_rates']['gt'] = float(regex.match(self.input_text).group(10)) # gt
            return True
        else:
            return False


    def _parse_raxml_aa(self):
        """
        Parse RAxML info file - amino acid.
        """

        # NH: there's really nothing wrong with the way you've
        # implemented any of these private methods - the only thing
        # that is significantly less efficient than it could be is
        # that the groups can be stored in an object - it isn't
        # necessary to re-run the regular expression for each group,
        # eg:
        # mobj = regex.match(self.input_text)
        # if mobj:
        #     program = mobj.groups()[1]
        #
        #
        # Also, note that regex.search (as opposed to regex.match)
        # allows matches starting in the middle of the string.
        #
        # however, I probably would have broken this up into multiple
        # regular expressions to make it easier to read and less
        # fragile.

        regex = re.compile('.*(RAxML version .*?) release.*DataType: (\w+).*Substitution Matrix: (\w+).*alpha\[0\]: (\d+\.\d+)', re.M|re.DOTALL)

        mobj = regex.match(self.input_text)
        if mobj:
            groups = mobj.groups()
            self.stats_values['program'] = groups[0] # program
            self.stats_values['datatype'] = groups[1] # datatype
            self.stats_values['subs_model'] = groups[2] # subs_model
            self.stats_values['ras_model'] = 'gamma' # gamma is present
            self.stats_values['gamma']['alpha'] = float(groups[3]) # alpha
            self.stats_values['gamma']['n_cats'] = 4 # n_cats - default is 4 for RaxML
            return True
        else:
            return False

    def _parse_phyml_dna(self):
        """
        Parse PhyML output file - DNA
        """

        regex = re.compile('.*---\s+(PhyML .*?)\s+ ---.*Model of nucleotides substitution:\s+(\w+).*Number of categories:\s+(\d+).*Gamma shape parameter:\s+(\d+\.\d+).*A <-> C\s+(\d+\.\d+).*A <-> G\s+(\d+\.\d+).*A <-> T\s+(\d+\.\d+).*C <-> G\s+(\d+\.\d+).*C <-> T\s+(\d+\.\d+).*G <-> T\s+(\d+\.\d+)',
                           re.M|re.DOTALL)
        #re_datatype = re.compile('.*Nucleotides frequencies.*', re.M|re.DOTALL)

        if (regex.match(self.input_text)):
            self.stats_values['program'] = regex.match(self.input_text).group(1) # program
            self.stats_values['subs_model'] = regex.match(self.input_text).group(2) # subs_model
            self.stats_values['ras_model'] = 'gamma' # gamma is present
            self.stats_values['gamma']['n_cats'] = int(regex.match(self.input_text).group(3)) # n_cats - default is 4 for RaxML
            self.stats_values['gamma']['alpha'] = float(regex.match(self.input_text).group(4)) # alpha
            self.stats_values['subs_rates']['ac'] = float(regex.match(self.input_text).group(5)) # ac
            self.stats_values['subs_rates']['ag'] = float(regex.match(self.input_text).group(6)) # ag
            self.stats_values['subs_rates']['at'] = float(regex.match(self.input_text).group(7)) # at
            self.stats_values['subs_rates']['cg'] = float(regex.match(self.input_text).group(8)) # cg
            self.stats_values['subs_rates']['ct'] = float(regex.match(self.input_text).group(9)) # ct
            self.stats_values['subs_rates']['gt'] = float(regex.match(self.input_text).group(10)) # gt
            self.stats_values['datatype'] = 'DNA' # datatype
            return True
        else:
            return False


    def _parse_phyml_aa(self):
        """
        Parse PhyML output file - AA
        """

        regex = re.compile('.*---\s+(PhyML .*?)\s+ ---.*Model of amino acids substitution:\s+(\w+).*Discrete gamma model:\s+(\w+)',
                           re.M|re.DOTALL)

        if (regex.match(self.input_text)):
            self.stats_values['program'] = regex.match(self.input_text).group(1) # program
            self.stats_values['subs_model'] = regex.match(self.input_text).group(2) # subs_model
            gamma_model = regex.match(self.input_text).group(3)
            if gamma_model.lower() == 'yes':
                self.stats_values['ras_model'] = 'gamma'
            else:
                self.stats_values['ras_model'] = None # gamma is not present
            self.stats_values['datatype'] = 'AA' # datatype
            return True
        else:
            return False


def build_parser(parser):
    parser.add_argument("-a", "--author",
        action="store", dest="author",
        help='Person who created the reference package', metavar='NAME')

    parser.add_argument("-c", "--clobber",
        action="store_true", dest="clobber", default = False,
        help= 'Delete an existing reference package.')
    
    parser.add_argument("-d", "--description",
        action="store", dest="description",
        help='An arbitrary description field', metavar='TEXT')

    parser.add_argument("-f", "--aln-fasta",
        action="store", dest="aln_fasta",
        help='Multiple alignment in fasta format', metavar='FILE')

    parser.add_argument("-i", "--seq-info",
        action="store", dest="seq_info",
        help='CSV format file describing the aligned reference ' + \
             'sequences, minimally containing the fields "seqname" ' + \
             'and "tax_id"', metavar='FILE')

    parser.add_argument("-l", "--locus",
        action="store", dest="locus", required=True,
        help='The locus described by the reference package', metavar='LOCUS')

    parser.add_argument("-m", "--mask",
        action="store", dest="mask",
        help='Text file containing a mask', metavar='FILE')

    parser.add_argument("-p", "--profile",
        action="store", dest="profile",
        help='Alignment profile', metavar='FILE')

    parser.add_argument('-P', '--package-name',
        action='store', dest='package_name',
        default='./taxtastic.refpkg', metavar='PATH',
        help='Name of output directory [default %(default)s]')

    parser.add_argument("-r", "--package-version",
        action="store", dest="package_version",
        help='Release version for the reference package', metavar='VERSION')

    parser.add_argument("-s", "--tree-stats",
        action="store", dest="tree_stats",
        help='File containing tree statistics (for example ' + \
             'RAxML_info.whatever")', metavar='FILE')

    parser.add_argument("-S", "--aln-sto",
        action="store", dest="aln_sto",
        help='Multiple alignment in Stockholm format', metavar='FILE')

    parser.add_argument("-t", "--tree-file",
        action="store", dest="tree_file",
        help='Phylogenetic tree in newick format',
        metavar='FILE')

    parser.add_argument("-T", "--taxonomy",
        action="store", dest="taxonomy",
        help='CSV format file defining the taxonomy. Fields include ' + \
             '"tax_id","parent_id","rank","tax_name" followed by a column ' + \
             'defining tax_id at each rank starting with root', metavar='FILE')

        
def action(args):
    """
    Create the reference package (a directory with a manifest named
    `manifest_name`).

     * args - output of argparse.ArgumentParser.parse_args (or
       presumably any other object with the required attributes:
       [TODO: list required attributes here])
     * MANIFEST_NAME - name of the JSON-format manifest file. Uses
       taxonomy.package.MANIFEST_NAME by default.
     * PACKAGE_CONTENTS - A dict defining sections and contents of
       each. Uses taxonomy.package.PACKAGE_CONTENTS by default.
     * PHYLO_MODEL_FILE - Names the JSON format file containing the
       PHYLO MODEL DATA; uses taxonomy.package.PHYLO_MODEL_FILE by default.
    """

    pkg_dir = args.package_name

    if args.clobber:
        utils.rmdir(pkg_dir)

    os.mkdir(pkg_dir)

    manifest = os.path.join(pkg_dir, MANIFEST_NAME)
    optdict = defaultdict(dict)

    # Add fields and values for the metadata section.
    optdict['metadata']['create_date'] = time.strftime('%Y-%m-%d %H:%M:%S')
    # locus is required
    optdict['metadata']['locus'] = args.locus
    # format_version is a module-level constant
    optdict['metadata']['format_version'] = FORMAT_VERSION
    if args.description:
        optdict['metadata']['description'] = args.description
    if args.author:
        optdict['metadata']['author'] = args.author
    if args.package_version:
        optdict['metadata']['package_version'] = args.package_version

    # phylo_model_file is part of the package, but not a command-line
    # argument; write out the phylo model file in JSON format, but
    # only if tree_stats was specified as an argument.
    if args.tree_stats:
        phylo_model_pth = os.path.join(pkg_dir, PHYLO_MODEL_FILE)

        parser = StatsParser(args.tree_stats)
        success = parser.parse_stats_data()

        if not success:
            raise ConfigError("Unable to create %s from %s." % (phylo_model_pth, args.tree_stats))

        parser.write_stats_json(phylo_model_pth)

        optdict['files']['phylo_model_file'] = PHYLO_MODEL_FILE
        optdict['md5']['phylo_model_file'] = \
            hashlib.md5(open(phylo_model_pth).read()).hexdigest()

    # copy all provided files into the package directory
    for fname in PACKAGE_CONTENTS['files']:
        if fname == 'phylo_model_file':
            continue

        pth = getattr(args, fname)
        if pth:
            shutil.copy(pth, pkg_dir)
            optdict['files'][fname] = os.path.split(pth)[1]
            optdict['md5'][fname] = hashlib.md5(open(pth).read()).hexdigest()

    package.write_config(fname=manifest, optdict=optdict)


