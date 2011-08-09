"""Reroots a reference package"""
import subprocess
import tempfile
import logging
import shutil
import os

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-p', '--pretend',
        action='store_true', default=False,
        help="don't save the rerooted tree; just attempt the rerooting.")

def action(args):
    refpkg_name, = args.refpkg
    rp = refpkg.Refpkg(refpkg_name)
    fd, name = tempfile.mkstemp(dir=refpkg_name)
    os.close(fd)
    try:
        subprocess.check_call(['rppr', 'reroot', '-c', refpkg_name, '-o', name])
    except:
        os.unlink(name)
        raise

    if args.pretend:
        os.unlink(name)
        return

    shutil.copystat(name, rp.resource_path('tree'))
    os.rename(name, rp.resource_path('tree'))
    rp.rehash('tree')
    rp.save()
