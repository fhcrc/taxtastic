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


def __version__():
    import pkg_resources
    def _safe_int(s):
        try:
            return int(s)
        except:
            return s
    version = pkg_resources.require("taxtastic")[0].version
    version_info = tuple(_safe_int(i) for i in  version.split('.'))
    return version, version_info

__version__, __version_info__ = __version__()
