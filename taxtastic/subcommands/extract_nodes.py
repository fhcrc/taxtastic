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
"""Extract nodes from a given source in yaml format

"""

import sys
import logging
import sqlalchemy
from itertools import groupby
from operator import itemgetter
# from collections import OrderedDict

import yaml
from fastalite import Opener

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import add_database_args

log = logging.getLogger(__name__)


def clean_dict(keys, vals):
    # outdict = OrderedDict()
    outdict = {}
    for k, v in zip(keys, vals):
        if k == 'source_id' or v is None:
            continue
        outdict[k] = bool(v) if k.startswith('is_') else v

    return outdict


def build_parser(parser):
    parser = add_database_args(parser)
    parser.add_argument('source_name',
                        help='name of source identifying names and nodes to extract')
    parser.add_argument('-o', '--outfile', type=Opener('w'), default=sys.stdout)


def action(args):
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 2)
    tax = Taxonomy(engine, schema=args.schema)

    with engine.connect() as con:
        # TODO: need to order nodes so that parents are always created first

        cmd = """
        select nodes.*, source.name as source_name
        from {nodes}
        join {source} on nodes.source_id = source.id
        where source.name = {x}
        """.format(x=tax.placeholder, nodes=tax.nodes, source=tax.source)

        result = con.execute(cmd, (args.source_name,))
        keys = list(result.keys())
        nodes = [clean_dict(keys, vals) for vals in result.fetchall()]

        # get the complete lineage for each node, and provide an
        # ordering for all nodes so that children may be placed after
        # parents.
        tax_ids = list(map(itemgetter('tax_id'), nodes))
        lineages = tax._get_lineage_table(tax_ids)
        ordering = {}
        for i, lineage in enumerate(lineages):
            tax_id = lineage[1]
            if tax_id not in ordering:
                ordering[tax_id] = i

        nodes = sorted(nodes, key=lambda n: ordering[n['tax_id']])

        cmd = """
        select names.*, source.name as source_name
        from {names}
        join {source} on names.source_id = source.id
        where source.name = {x}
        """.format(x=tax.placeholder, names=tax.names, source=tax.source)

        result = con.execute(cmd, (args.source_name,))
        keys = list(result.keys())
        names = [clean_dict(keys, vals) for vals in result.fetchall()]
        namedict = {key: list(grp)
                    for key, grp in groupby(names, itemgetter('tax_id'))}

        for node in nodes:
            node['type'] = 'node'
            tax_id = node['tax_id']
            if tax_id in namedict:
                node['names'] = namedict.pop(tax_id)

        yaml.safe_dump_all(nodes, args.outfile, default_flow_style=False,
                           explicit_start=True, indent=2)

        # prepare remaining names
        remaining_names = []
        for tax_id, names in list(namedict.items()):
            for name in names:
                del name['tax_id']

            remaining_names.append({
                'tax_id': tax_id,
                'type': 'name',
                'names': names
            })

        yaml.safe_dump_all(remaining_names, args.outfile, default_flow_style=False,
                           explicit_start=True, indent=2)
