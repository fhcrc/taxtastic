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
import glob
from os.path import splitext, split, join, dirname

commands = [
    'add_nodes',
    'check',
    'info',
    'create',
    'new_database',
    'reroot',
    'update',
    'taxids',
    'update_taxids',
    'taxtable',
    'strip',
    'rollback',
    'rollforward',
    'rp',
    'refpkg_intersection',
    'merge'
]


def itermodules(subcommands_path=None, root=__name__):

    if subcommands_path is None:
        subcommands_path = dirname(__file__)

    modules = sorted(glob.glob(join(subcommands_path, '*.py')))
    # excluded = set(['lonelynodes'])
    excluded = set()

    commands = [x for x in [splitext(split(p)[1])[0]
                            for p in modules]
                if not x.startswith('_') and x not in excluded]

    for command in commands:
        yield command, __import__('%s.%s' % (root, command), fromlist=[command])


def close_all_files(args):
    for name, obj in args.__dict__.items():
        if obj and hasattr(obj, 'close'):
            obj.close()
