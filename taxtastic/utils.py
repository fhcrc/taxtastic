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
import datetime
import logging
import os
import re
import subprocess

log = logging

try:
    # download: http://pypi.python.org/pypi/xlrd
    # docs: http://www.lexicon.net/sjmachin/xlrd.html
    import xlrd
except ImportError:
    xlrd = None

if xlrd:
    def _cellval(cell_obj, datemode):
        if cell_obj.ctype == xlrd.XL_CELL_DATE:
            timetup = xlrd.xldate_as_tuple(cell_obj.value, datemode)
            val = datetime.datetime(*timetup)
        elif cell_obj.ctype == xlrd.XL_CELL_TEXT:
            # coerce to pain text from unicode
            val = str(cell_obj.value).strip()
        else:
            val = cell_obj.value

        return val

    def read_spreadsheet(filename, fmts=None):
        """
        Read excel spreadsheet, performing type coersion as specified
        in fmts (dict keyed by column name returning either a
        formatting string or a function such as str, int, float,
        etc). Returns (list headers, iter rows)
        """

        w = xlrd.open_workbook(filename)
        datemode = w.datemode
        s = w.sheet_by_index(0)
        rows = ([_cellval(c, datemode) for c in s.row(i)] for i in xrange(s.nrows))

        firstrow = rows.next()
        headers = [str('_'.join(x.split())) for x in firstrow]

        lines = []
        for row in rows:
            # valid rows have at least one value
            if not any([bool(cell) for cell in row]):
                continue

            d = dict(zip(headers, row))

            if fmts:
                for colname in fmts.keys():
                    if hasattr(fmts[colname], '__call__'):
                        formatter = fmts[colname]
                    else:
                        formatter = lambda val: fmts[colname] % val

                    try:
                        d[colname] = formatter(d[colname])
                    except (TypeError, ValueError, AttributeError):
                        pass

            lines.append(d)

        return headers, iter(lines)

def get_new_nodes(fname):
    """
    Return an iterator of dicts given either an .xls spreadsheet or
    .csv-format file.
    """

    if fname.endswith('.xls'):
        if not xlrd:
            raise AttributeError('xlrd not installed: cannot parse .xls files.')

        fmts = {'tax_id':'%i', 'parent_id':'%i'}
        headers, rows = read_spreadsheet(fname, fmts)
    elif fname.endswith('.csv'):
        with open(fname, 'rU') as infile:
            reader = list(csv.DictReader(infile))
            rows = (d for d in reader if d['tax_id'])
    else:
        raise ValueError('Error: %s must be in .csv or .xls format')


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

    with open(fname) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                yield line.split('#', 1)[0].strip()


def try_set_fields(d, regex, text, hook=lambda x: x):
    v = re.search(regex, text, re.MULTILINE)
    if v:
        d.update(dict([(key,hook(val)) for key,val
                       in v.groupdict().iteritems()]))
    return d

class InvalidLogError(ValueError):
    pass

def parse_raxml(handle):
    """Parse RAxML's summary output.

    *handle* should be an open file handle containing the RAxML
    output.  It is parsed and a dictionary returned.
    """
    s = ''.join(handle.readlines())
    result = {}
    try_set_fields(result, r'(?P<program>RAxML version [0-9.]+)', s)
    try_set_fields(result, r'(?P<datatype>DNA|RNA|AA)', s)
    result['empirical_frequencies'] = not(result['datatype'] == 'AA') or \
        re.search('Empirical Base Frequencies', s) != None
    try_set_fields(result, r'Substitution Matrix: (?P<subs_model>\w+)', s)
    rates = {}
    try_set_fields(rates,
                   (r"rates\[0\] ac ag at cg ct gt: "
                    r"(?P<ac>[0-9.]+) (?P<ag>[0-9.]+) (?P<at>[0-9.]+) "
                    r"(?P<cg>[0-9.]+) (?P<ct>[0-9.]+) (?P<gt>[0-9.]+)"), s, hook=float)
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


JTT_MODEL = 'ML Model: Jones-Taylor-Thorton, CAT approximation with 20 rate categories'
WAG_MODEL = 'ML Model: Whelan-And-Goldman, CAT approximation with 20 rate categories'

def parse_fasttree(fobj):
    data = {
        'empirical_frequencies': True,
        'datatype': 'DNA',
        'subs_model': 'GTR',
        'ras_model': 'Price-CAT',
        'Price-CAT': {},
    }
    for line in fobj:
        if not line: continue
        splut = line.split()
        if splut[0] == 'FastTree':
            data['program'] = line.strip()
        elif splut[0] == 'Rates':
            data['Price-CAT']['Rates'] = map(float, splut[1:])
        elif splut[0] == 'SiteCategories':
            data['Price-CAT']['SiteCategories'] = map(int, splut[1:])
        elif splut[0] == 'NCategories':
            data['Price-CAT']['n_cats'] = int(splut[1])
        elif splut[0] == 'GTRRates':
            data['subs_rates'] = dict(
                zip(['ac', 'ag', 'at', 'cg', 'ct', 'gt'],
                    map(float, splut[1:])))
        elif line.strip() == JTT_MODEL:
            data['subs_model'] = 'JTT'
            data['datatype'] = 'AA'
        elif line.strip() == WAG_MODEL:
            data['subs_model'] = 'WAG'
            data['datatype'] = 'AA'

    # Sanity check
    if data['subs_model'] == 'GTR' and 'subs_rates' not in data:
        raise InvalidLogError("GTR model, but no substitution rates found!")

    return data

def parse_phyml(fobj):
    s = ''.join(fobj)
    result = {'gamma': {}}
    try_set_fields(result, r'---\s*(?P<program>PhyML.*?)\s*---', s)
    try_set_fields(result['gamma'], r'Number of categories:\s+(?P<n_cats>\d+)',
                   s, hook=int)
    try_set_fields(result['gamma'],
            r'Gamma shape parameter:\s+(?P<alpha>\d+\.\d+)', s, hook=float)
    result['ras_model'] = 'gamma'
    if 'nucleotides' in s:
        result['datatype'] = 'DNA'
        try_set_fields(result,
                       r'Model of nucleotides substitution:\s+(?P<subs_model>\w+)',
                       s)
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
        try_set_fields(result,
                       r'Model of amino acids substitution:\s+(?P<subs_model>\w+)',
                       s)
    else:
        raise ValueError('Could not determine if alignment is AA or DNA')

    return result

def has_rppr(rppr_name='rppr'):
    """
    Check for rppr binary in path
    """
    with open(os.devnull) as dn:
        try:
            subprocess.check_call([rppr_name], stdout=dn, stderr=dn)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                return False
            else:
                raise
        except subprocess.CalledProcessError as e:
            # rppr returns non-zero exit status with no arguments
            pass
    return True
