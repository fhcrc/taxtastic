import datetime
import logging
import shutil
import os
import csv
from os import path

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
                    except (TypeError, ValueError, AttributeError), msg:
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

    return rows


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

def mkdir(dirpath, clobber = False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber:
        rmdir(dirpath)

    try:
        os.mkdir(dirpath)
    except OSError, msg:
        log.debug(msg)

    if not path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath

def rmdir(dirpath):
    """
    Remove a (potentially missing) directory without errors. Raise
    OSError if directory can't be removed.
    """

    try:
        shutil.rmtree(dirpath)
    except OSError, msg:
        log.debug(msg)

    if path.exists(dirpath):
        raise OSError('Failed to remove %s' % dirpath)


