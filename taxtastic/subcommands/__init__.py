commands = [
    'add_nodes',
    'badgraph',
    'check',
    'convexify',
    'create',
    'new_database',
    'reroot',
    'taxids',
    'taxtable',
    ]

def itermodules(root=__name__):
    for command in commands:
        yield command, __import__('%s.%s' % (root, command), fromlist=[command])
