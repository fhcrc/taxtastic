"""Not yet implemented"""

def build_parser(parser):
    parser.add_argument('-P', '--package-name',
        action='store', dest='package_name',
        default='./taxtastic.refpkg', metavar='PATH',
        help='Name of output directory [default %(default)s]')

def action(args):
    return 1
