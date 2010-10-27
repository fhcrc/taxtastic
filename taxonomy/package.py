import ConfigParser
import logging
import os
import time
import shutil
import hashlib

log = logging

PACKAGE_VERSION = 0

manifest_name = 'CONTENTS.json'

package_contents = {
    'metadata':['create_date','author','description','package_version'],
    'files':['tree_file','tree_stats','aln_fasta','aln_sto',
             'profile', 'seq_info','taxonomy','mask'],
    'md5':[]
    }

def write_config(fname, optdict, sections):
    """
    * fname - name of config file
    * optdict - a dict with keys d[(section, opt)] = val
    * sections - a list of tuples, eg
      (('section1',['opt1','opt2',...]),('section2',['opt3','opt4',...]))
      where 'opt*' are attributes of 'options'
    """

    config = ConfigParser.SafeConfigParser()

    for section_name, opts in sections.items():
        config.add_section(section_name)
        for opt in opts:
            val = optdict.get((section_name,opt),'')
            config.set(section_name, opt, str(val))

    cfile = open(fname,'w+')
    config.write(cfile)
    log.info('writing %s' % fname)

    cfile.seek(0)
    # log.info('\n' + cfile.read())
    print cfile.read()
    cfile.close()

def create(pkg_dir, options, manifest_name=manifest_name, package_contents=package_contents):

    os.mkdir(pkg_dir)
    manifest = os.path.join(pkg_dir, manifest_name)

    optdict = {}
    optdict[('metadata','create_date')] = time.strftime('%Y-%m-%d %H:%M:%S')

    # copy files into the package directory
    for fname in package_contents['files']:
        pth = getattr(options, fname)

        if pth:
            shutil.copy(pth, pkg_dir)
            optdict[('files',fname)] = os.path.split(pth)[1]
            optdict[('md5',fname)] = hashlib.md5(open(pth).read()).hexdigest()
            package_contents['md5'].append(fname)

    write_config(fname=manifest, optdict=optdict, sections=package_contents)





# StatsParser class - parse tree stats output files generated 
# by RAxML (2 types) and PhyML.
#
# Supported file types include:
#     raxml_condensed - subs_rates on a single line
#     raxml_re_estimated - subs_rates 1 per line
#     raxml_aa - amino acid
#     phyml_dna
#     phyml_aa

import os
import re
import json
from collections import defaultdict

class StatsParser(object):
    """Tree Statistics Output Parser Class"""
    
    file_name = '' # Name of file to parse.
    file_type = 'unknown'  # Type of file examined, defaults to unknown. 
    input_text = '' # Full text of input file.
    # Property of the class that gets populated when matches are found 
    # within an output file.
    stats_values = defaultdict(dict)

### Public methods ###    
    
    # Return contents of stats_values, in the format of a dictionary.
    def get_stats_values(self):
        return self.stats_values
    
    # Return contents of stats_values in JSON format.
    def get_stats_json(self):
        return json.dumps(self.stats_values, indent=2)
    
    def write_stats_json(self, out_file):
        try:
            json = self.get_stats_json()
            out = open(out_file, 'w')
            out.write(json)
            out.close()
        except IOError as (errno, strerror):
            raise Exception, "I/O error({0}): {1}".format(errno, strerror)    
    
    # Return contents of file_type.
    def get_file_type(self):
        return self.file_type
     
    # Extract data from stats output files.
    def parse_stats_data(self):    
        # Read input file into a string for multiline matching.
        match_found = False
        try:
            input_file = open(self.file_name, 'r')
            self.input_text = input_file.read()
            input_file.close()
        except IOError as (errno, strerror):
            raise Exception, "I/O error({0}): {1}".format(errno, strerror)
        if (self._parse_raxml_condensed()):
            self.file_type = 'raxml_condensed'
            match_found = True
        elif (self._parse_raxml_re_estimated()):
            self.file_type = 'raxml_re_estimated'
            match_found = True            
        elif (self._parse_raxml_aa()):
            self.file_type = 'raxml_aa'
            match_found = True            
        elif (self._parse_phyml_dna()):
            self.file_type = 'phyml_dna'
            match_found = True       
        elif (self._parse_phyml_aa()):
            self.file_type = 'phyml_aa'
            match_found = True   
        return match_found    
    
    # Determine the md5sum of the output file.
    def get_md5sum(self):
        pass
    
### Private methods ###    
    
    # Parse RAxML info file - condensed.    
    def _parse_raxml_condensed(self):
        try:
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
        except:
            raise Exception, "Encountered a problem within  _parse_raxml_condensed."

    # Parse RAxML info file - re-estimated.
    def _parse_raxml_re_estimated(self):
        regex = re.compile('.*(RAxML version .*?) release.*DataType: (\w+).*Substitution Matrix: (\w+).*alpha: (\d+\.\d+).*rate A <-> C: (\d+\.\d+).*rate A <-> G: (\d+\.\d+).*rate A <-> T: (\d+\.\d+).*rate C <-> G: (\d+\.\d+).*rate C <-> T: (\d+\.\d+).*rate G <-> T: (\d+\.\d+)', 
                           re.M|re.DOTALL)
        try:
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
        except:
            raise Exception, "Encountered a problem within _parse_raxml_re_estimated."
     

    # Parse RAxML info file - amino acid.
    def _parse_raxml_aa(self):
        regex = re.compile('.*(RAxML version .*?) release.*DataType: (\w+).*Substitution Matrix: (\w+).*alpha\[0\]: (\d+\.\d+)',
                           re.M|re.DOTALL)
        try:
            if (regex.match(self.input_text)):
                self.stats_values['program'] = regex.match(self.input_text).group(1) # program
                self.stats_values['datatype'] = regex.match(self.input_text).group(2) # datatype
                self.stats_values['subs_model'] = regex.match(self.input_text).group(3) # subs_model
                self.stats_values['ras_model'] = 'gamma' # gamma is present
                self.stats_values['gamma']['alpha'] = float(regex.match(self.input_text).group(4)) # alpha
                self.stats_values['gamma']['n_cats'] = 4 # n_cats - default is 4 for RaxML
                return True       
            else:
                return False   
        except:
            raise Exception, "Encountered a problem within _parse_raxml_aa."
     




    # Parse PhyML output file - DNA
    def _parse_phyml_dna(self):
        regex = re.compile('.*---\s+(PhyML .*?)\s+ ---.*Model of nucleotides substitution:\s+(\w+).*Number of categories:\s+(\d+).*Gamma shape parameter:\s+(\d+\.\d+).*A <-> C\s+(\d+\.\d+).*A <-> G\s+(\d+\.\d+).*A <-> T\s+(\d+\.\d+).*C <-> G\s+(\d+\.\d+).*C <-> T\s+(\d+\.\d+).*G <-> T\s+(\d+\.\d+)', 
                           re.M|re.DOTALL)
        re_datatype = re.compile('.*Nucleotides frequencies.*', re.M|re.DOTALL)
        try:
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
        except:
            raise Exception, "Encountered a problem within _parse_phyml_dna."

    # Parse PhyML output file - AA
    def _parse_phyml_aa(self):
        regex = re.compile('.*---\s+(PhyML .*?)\s+ ---.*Model of amino acids substitution:\s+(\w+).*', 
                           re.M|re.DOTALL)

        try:
            if (regex.match(self.input_text)):
                self.stats_values['program'] = regex.match(self.input_text).group(1) # program
                self.stats_values['subs_model'] = regex.match(self.input_text).group(2) # subs_model
                self.stats_values['ras_model'] = None # gamma is not present
                self.stats_values['datatype'] = 'AA' # datatype   
                return True       
            else:
                return False   
        except:
            raise Exception, "Encountered a problem within _parse_phyml_aa."

            
### Constructor ###
            
    def __init__(self, file_name):
        self.file_name = file_name
        # Make sure the file passed in both exists and is readable.
        try:            
            if (not os.access(self.file_name, os.R_OK)):
                raise Exception()
        except:
            raise Exception, 'File is not readable: ' + file_name
