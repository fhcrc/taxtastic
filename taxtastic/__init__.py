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


__version__ = '0.5.0'
__version_info__ = (0, 5, 0)

if '+' in __version__ or 'pre' in __version__:
    # Try to append the commit hash to the version
    try:
        import subprocess
        p = subprocess.Popen(['git', 'log', '--pretty=format:%h', '-n', '1'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if out:
            __version__ += '-' + out.strip()
    except Exception:
        pass
