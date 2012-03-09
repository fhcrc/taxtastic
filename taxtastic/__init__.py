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

from os import path

def _safeint(s):
    try:
        return int(s)
    except ValueError:
        return s

_shafile = path.join(path.dirname(__file__), 'data', 'sha')
try:
    with open(_shafile) as f:
        sha = f.read().strip()
except Exception, e:
    sha = ''

__version__ = "0.3.2" + ('.' + sha if sha else '')
__version_info__ = tuple([_safeint(num) for num in __version__.split('.')])



