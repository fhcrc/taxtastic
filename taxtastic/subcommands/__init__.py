commands = [
    'add_nodes',
    'badgraph',
    'check',
    'create',
    'new_database',
    'reroot',
    'update',
    'taxids',
    'taxtable',
    'strip',
    'rollback',
    'rollforward',
    ]

def itermodules(root=__name__):
    for command in commands:
        yield command, __import__('%s.%s' % (root, command), fromlist=[command])
