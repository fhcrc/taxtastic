import datetime
import logging

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
            val = str(cell_obj.value)
        else:
            val = cell_obj.value

        return val

    def read_spreadsheet(filename, fmts=None):
        """
        Read excel spreadsheet, performing type coersion as specified in
        fmts (dict keyed by column name returning either a formatting
        string or a function such as str, int, float, etc).
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

        return headers, lines

def get_new_nodes(fname):
    if fname.endswith('.xls'):
        if xlrd:
            headers, rows = read_spreadsheet(
                fname,
                fmts={'tax_id':'%i', 'parent_id':'%i' }
                )
        else:
            raise AttributeError('xlrd not installed: cannot parse .xls files.')
    elif fname.endswith('.csv'):
        reader = csv.DictReader(fname)
        rows = (d for d in reader if d['tax_id'])
    else:
        raise ValueError('Error: %s must be in .csv or .xls format')

    return rows


