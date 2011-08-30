"""Reroots a reference package"""
import subprocess
import tempfile
import logging
import shutil
import os

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('--rppr', action='store', default=None,
                        help="specify the rppr binary to call to perform the rerooting")
    parser.add_argument('-p', '--pretend',
                        action='store_true', default=False,
                        help="don't save the rerooted tree; just attempt the rerooting.")

def action(args):
    r = refpkg.Refpkg(args.refpkg)
    r.reroot(rppr=args.rppr, pretend=args.pretend)
