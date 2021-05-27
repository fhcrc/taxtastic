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
import csv
import errno
import logging
import os
import re
import subprocess
import string
import random
from collections import OrderedDict

from six.moves import configparser

log = logging


def get_new_nodes(fname):
    """
    Return an iterator of dicts given a .csv-format file.
    """

    with open(fname, 'rU') as infile:
        infile = (line for line in infile if not line.startswith('#'))
        reader = list(csv.DictReader(infile))
        rows = (d for d in reader if d['tax_id'])

    # for now, children are provided as a semicolon-delimited list
    # within a cell (yes, yuck). We need to convert thit into a list
    # if present.
    for d in rows:
        if 'children' in d:
            if d['children']:
                d['children'] = [x.strip() for x in d['children'].split(';')]
            else:
                del d['children']
        yield d


def getlines(fname):
    """
    Returns iterator of whitespace-stripped lines in file, omitting
    blank lines, lines beginning with '#', and line contents following
    the first '#' character.
    """

    with open(fname, 'rU') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                yield line.split('#', 1)[0].strip()


def try_set_fields(d, regex, text, hook=lambda x: x):
    v = re.search(regex, text, re.MULTILINE)
    if v:
        d.update(dict([(key, hook(val)) for key, val
                       in v.groupdict().items()]))
    return d


class InvalidLogError(ValueError):
    pass

def parse_raxmlng(handle):
    """Parse RAxMLng's summary output.

    *handle* should be an open file handle containing the RAxMLng
    log.  It is parsed and a dictionary returned.
    """
    s = handle.read()
    result = {}
    try_set_fields(result, r'(?P<program>RAxML-NG v. [\d\.]+)', s)
    # Less ideal, but for now force DNA and GTR for now
    result['datatype'] = 'DNA'
    result["subs_model"] =  "GTR"
    try_set_fields(result, r'\\nModel: (?P<subs_model>[\w\+]+)\\n', s)
    result['empirical_frequencies'] = re.search(r'Base frequencies \(ML\)', s) is not None
    rates = {}
    try_set_fields(rates, r'Substitution rates \(ML\): (?P<ac>\d+\.\d+) (?P<ag>\d+\.\d+) (?P<at>\d+\.\d+) (?P<cg>\d+\.\d+) (?P<ct>\d+\.\d+) (?P<gt>\d+\.\d+)', s, hook=float)
    if len(rates) > 0:
        result['subs_rates'] = rates
    gamma = {}
    try_set_fields(gamma, r'Rate heterogeneity: GAMMA \((?P<n_cats>\d+) cats, mean\),  alpha: (?P<alpha>\d+\.\d+)', s, hook=float)
    try:
        gamma['n_cats'] = int(gamma['n_cats'])
    except:
        pass
    result['gamma'] = gamma
    result['ras_model'] = 'gamma'

    return result        

# https://github.com/amkozlov/raxml-ng/wiki/Input-data#single-mode
DNA = ['JC', 'K80', 'F81', 'HKY', 'TN93ef', 'TN93', 'K81', 'K81uf', 'TPM2',
       'TPM2uf', 'TPM3', 'TPM3uf', 'TIM1', 'TIM1uf', 'TIM2', 'TIM2uf',
       'TIM3', 'TIM3uf', 'TVMef', 'TVM', 'SYM', 'GTR']
PROT = ['Blosum62', 'cpREV', 'Dayhoff', 'DCMut', 'DEN', 'FLU', 'HIVb',
        'HIVw', 'JTT', 'JTT-DCMut', 'LG', 'mtART', 'mtMAM', 'mtREV',
        'mtZOA', 'PMB', 'rtREV', 'stmtREV', 'VT', 'WAG', 'LG4M', 'LG4X',
        'PROTGTR']


def parse_raxmlng(handle):
    """Parse RAxMLng's summary output.
    *handle* should be an open file handle containing the RAxMLng
    log.  It is parsed and a dictionary returned.
    """
    s = handle.read()
    result = {}
    try_set_fields(result, r'(?P<program>RAxML-NG v. [\d\.]+)', s)
    try_set_fields(result, r'Model: (?P<subs_model>\w+)', s)
    try_set_fields(result, r'raxml-ng.*?--data-type (?P<datatype>\w+)', s)
    if 'datatype' not in result:
        if result['subs_model'] in DNA:
            result['datatype'] = 'DNA'
        elif result['subs_model'] in PROT:
            result['datatype'] = 'RNA'  # AA specified in --data-type
        else:
            result['datatype'] = ''  # There's more but stop here for now
    result['empirical_frequencies'] = re.search(r'Base frequencies \(ML\)', s) is not None
    rates = {}
    try_set_fields(rates, r'Substitution rates \(ML\): (?P<ac>\d+\.\d+) (?P<ag>\d+\.\d+) (?P<at>\d+\.\d+) (?P<cg>\d+\.\d+) (?P<ct>\d+\.\d+) (?P<gt>\d+\.\d+)', s, hook=float)
    if len(rates) > 0:
        result['subs_rates'] = rates
    gamma = {}
    try_set_fields(gamma, r'Rate heterogeneity: GAMMA \((?P<n_cats>\d+) cats, mean\),  alpha: (?P<alpha>\d+\.\d+)', s, hook=float)
    if gamma:
        gamma['n_cats'] = int(gamma['n_cats'])
        result['gamma'] = gamma
        result['ras_model'] = 'gamma'
    return result


def parse_raxml(handle):
    """Parse RAxML's summary output.

    *handle* should be an open file handle containing the RAxML
    output.  It is parsed and a dictionary returned.
    """
    s = ''.join(handle.readlines())
    result = {}
    try_set_fields(result, r'(?P<program>RAxML version [0-9.]+)', s)
    try_set_fields(result, r'DataType: (?P<datatype>DNA|RNA|AA)', s)
    result['empirical_frequencies'] = (
        result['datatype'] != 'AA' or
        re.search('empirical base frequencies', s, re.IGNORECASE) is not None)
    try_set_fields(result, r'Substitution Matrix: (?P<subs_model>\w+)', s)
    rates = {}
    if result['datatype'] != 'AA':
        try_set_fields(rates,
                       (r"rates\[0\] ac ag at cg ct gt: "
                        r"(?P<ac>[0-9.]+) (?P<ag>[0-9.]+) (?P<at>[0-9.]+) "
                        r"(?P<cg>[0-9.]+) (?P<ct>[0-9.]+) (?P<gt>[0-9.]+)"),
                       s, hook=float)
        try_set_fields(rates, r'rate A <-> C: (?P<ac>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate A <-> G: (?P<ag>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate A <-> T: (?P<at>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate C <-> G: (?P<cg>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate C <-> T: (?P<ct>[0-9.]+)', s, hook=float)
        try_set_fields(rates, r'rate G <-> T: (?P<gt>[0-9.]+)', s, hook=float)
        if len(rates) > 0:
            result['subs_rates'] = rates
    result['gamma'] = {'n_cats': 4}
    try_set_fields(result['gamma'],
                   r"alpha[\[\]0-9]*: (?P<alpha>[0-9.]+)", s, hook=float)
    result['ras_model'] = 'gamma'
    return result


JTT_MODEL = ('ML Model: Jones-Taylor-Thorton, CAT '
             'approximation with 20 rate categories')
WAG_MODEL = ('ML Model: Whelan-And-Goldman, CAT '
             'approximation with 20 rate categories')
LG_MODEL = 'ML Model: Le-Gascuel 2008, CAT approximation with 20 rate categories'


def parse_fasttree(fobj):
    data = {
        'empirical_frequencies': True,
        'datatype': 'DNA',
        'subs_model': 'GTR',
        'ras_model': 'Price-CAT',
        'Price-CAT': {},
    }
    for line in fobj:
        if not line.strip():
            continue
        splut = line.split()
        if splut[0] == 'FastTree':
            data['program'] = line.strip()
        elif splut[0] == 'Rates':
            data['Price-CAT']['Rates'] = list(map(float, splut[1:]))
        elif splut[0] == 'SiteCategories':
            data['Price-CAT']['SiteCategories'] = list(map(int, splut[1:]))
        elif splut[0] == 'NCategories':
            data['Price-CAT']['n_cats'] = int(splut[1])
        elif splut[0] == 'GTRRates':
            data['subs_rates'] = dict(
                list(zip(['ac', 'ag', 'at', 'cg', 'ct', 'gt'],
                         list(map(float, splut[1:])))))
        elif line.strip() == JTT_MODEL:
            data['subs_model'] = 'JTT'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False
        elif line.strip() == WAG_MODEL:
            data['subs_model'] = 'WAG'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False
        elif line.strip() == LG_MODEL:
            data['subs_model'] = 'LG'
            data['datatype'] = 'AA'
            data['empirical_frequencies'] = False

    # Sanity check
    if data['subs_model'] == 'GTR' and 'subs_rates' not in data:
        raise InvalidLogError("GTR model, but no substitution rates found!")

    return data


def parse_phyml(fobj, frequency_type=None):
    if frequency_type not in (None, 'empirical', 'model'):
        raise ValueError("Unknown frequency_type: {0}".format(frequency_type))
    s = ''.join(fobj)
    result = {'gamma': {}}
    try_set_fields(result, r'---\s*(?P<program>PhyML.*?)\s*---', s)
    try_set_fields(
        result['gamma'],
        r'Number of categories:\s+(?P<n_cats>\d+)', s, hook=int)
    try_set_fields(
        result['gamma'],
        r'Gamma shape parameter:\s+(?P<alpha>\d+\.\d+)', s, hook=float)
    result['ras_model'] = 'gamma'
    if 'nucleotides' in s:
        result['datatype'] = 'DNA'
        try_set_fields(
            result,
            r'Model of nucleotides substitution:\s+(?P<subs_model>\w+)', s)
        rates = {}
        try_set_fields(rates, r'A <-> C\s+(?P<ac>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'A <-> G\s+(?P<ag>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'A <-> T\s+(?P<at>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'C <-> G\s+(?P<cg>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'C <-> T\s+(?P<ct>\d+\.\d+)', s, hook=float)
        try_set_fields(rates, r'G <-> T\s+(?P<gt>\d+\.\d+)', s, hook=float)
        if rates:
            result['subs_rates'] = rates

        # PhyML doesn't record whether empirical base frequencies were used, or
        # ML estimates were made.
        # Setting to empirical for now.
        result['empirical_frequencies'] = True
    elif 'amino acids' in s:
        result['datatype'] = 'AA'
        if not frequency_type:
            raise ValueError("frequency type required for PhyML AA models.")
        result['empirical_frequencies'] = frequency_type == 'empirical'
        try_set_fields(
            result,
            r'Model of amino acids substitution:\s+(?P<subs_model>\w+)', s)
    else:
        raise ValueError('Could not determine if alignment is AA or DNA')

    return result


def parse_stockholm(fobj):
    """Return a list of names from an Stockholm-format sequence alignment
    file. ``fobj`` is an open file or another object representing a
    sequence of lines.

    """

    names = OrderedDict()

    found_eof = False
    for line in fobj:
        line = line.strip()
        if line == '//':
            found_eof = True
        elif line.startswith('#') or not line.strip():
            continue
        else:
            name, __ = line.split(None, 1)
            names[name] = None

    if not found_eof:
        raise ValueError('Invalid Stockholm format: no file terminator')

    return list(names.keys())


def has_rppr(rppr_name='rppr'):
    """
    Check for rppr binary in path
    """
    with open(os.devnull) as dn:
        try:
            subprocess.check_call([rppr_name], stdout=dn, stderr=dn)
        except OSError as e:
            if e.errno == errno.ENOENT:
                return False
            else:
                raise
        except subprocess.CalledProcessError as e:
            # rppr returns non-zero exit status with no arguments
            pass
    return True


def add_database_args(parser):
    '''
    Add a standard set of database arguments for argparse
    '''
    parser.add_argument(
        'url',
        nargs='?',
        default='sqlite:///ncbi_taxonomy.db',
        type=sqlite_default(),
        help=('Database string URI or filename.  If no database scheme '
              'specified \"sqlite:///\" will be prepended. [%(default)s]'))
    db_parser = parser.add_argument_group(title='database options')

    # TODO: better description of what --schema does
    db_parser.add_argument(
        '--schema',
        help=('Name of SQL schema in database to query '
              '(if database flavor supports this).'))

    return parser


def sqlite_default():
    '''
    Prepend default scheme if none is specified. This helps provides backwards
    compatibility with old versions of taxtastic where sqlite was the automatic
    default database.
    '''
    def parse_url(url):
        # TODO: need separate option for a config file
        if url.endswith('.db') or url.endswith('.sqlite'):
            if not url.startswith('sqlite:///'):
                url = 'sqlite:///' + url
        elif url.endswith('.cfg') or url.endswith('.conf'):
            conf = configparser.SafeConfigParser(allow_no_value=True)
            conf.optionxform = str  # options are case-sensitive
            conf.read(url)
            url = conf.get('sqlalchemy', 'url')

        return url
    return parse_url


def random_name(length):
    return ''.join([random.choice(string.ascii_letters) for n in range(length)])
