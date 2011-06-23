#!/usr/bin/env python

"""
usage:

dev/set_version.py

Retrieve the curent version stored as an annotated git tag in the
format 'major.minor.release', and update bioscons/__init__.py to
reflect the new vaule. If you want to update the tag, do this:

git tag -a -m "increment version" major.minor.release

"""

from subprocess import Popen, PIPE
import fileinput
import sys

def _safeint(s):
    try:
        return int(s)
    except ValueError:
        return s

def shell(args):
    return Popen(args, stdout=PIPE).communicate()[0].strip()

def main():

    if '-h' in sys.argv:
        print __doc__
        sys.exit()

    dryrun = '-n' in sys.argv

    current_ver = shell(["git", "describe"])

    ver_info = [_safeint(x) for x in current_ver.split('.')]

    new_release = '.'.join('%s'%x for x in ver_info)
    new_version = '.'.join('%s'%x for x in ver_info[:2])

    target = 'taxtastic/__init__.py'

    print 'updating version in %s' % target
    for line in fileinput.input(target, inplace = not dryrun):
        if line.startswith('__version__'):
            line = '__version__ = "%s"\n' % new_release
        sys.stdout.write(line)
    fileinput.close()

if __name__ == '__main__':
    main()


