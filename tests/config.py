import os
from os import path, mkdir
import sys
import logging
import re
import unittest
import commands

log = logging

def funcname(idstr):
    return '.'.join(idstr.split('.')[1:])

from taxtastic import ncbi

# set verbosity of logging output
try:
    logflag = re.findall(r'-[vq]+\b', ' '.join(sys.argv[1:]))[0]
except IndexError:
    logflag = ''

logging.basicConfig(
    file = sys.stdout,
    format = '%(levelname)s %(module)s %(lineno)s %(message)s' \
        if logflag.startswith('-v') else '%(message)s',
    level = {'-q':logging.ERROR,
             '':logging.WARNING,
             '-v': logging.INFO,
             '-vv': logging.DEBUG}[logflag]
    )

# module data
datadir = path.join(path.dirname(__file__), '..', 'testfiles')
outputdir = path.join(path.dirname(__file__), '..', 'test_output')

if not os.path.isdir(outputdir):
    mkdir(outputdir)

class TestBase(unittest.TestCase):
    """
    Base class for unit tests
    """

    outputdir = outputdir

    def outdir(self):
        """
        Name an outputdir as outpudir/module.class.method
        """
        funcname = '.'.join(self.id().split('.')[-3:])
        return path.join(self.outputdir, funcname)

    def mkoutdir(self, clobber = True):
        """
        Create output directory (destructively if clobber is True)
        """

        outdir = self.outdir()
        mkdir(outdir, clobber)
        return outdir


class TestScriptBase(TestBase):

    """
    Base class for unit tests of scripts distributed with this
    package. Create a new test module like this::

      from config import TestScriptBase
      TestScriptBase.executable = './script_name.py'
      TestScriptBase.outputdir = 'path/to/somewhere'

      class TestSomething(TestScriptBase):

          def test01(self):
              self.cmd_ok('--help')

    """

    executable = None

    def __getitem__(self, i):
        """
        Enables string formatting, eg:
        print 'current function: %(funcname)s' % self
        """
        return getattr(self, i)

    def wrap_cmd(self, args = None, cmd = None):
        if cmd is None:
            cmd = self.executable
        input = (cmd + ' ' + args) % self
        log.warning('--> '+ input)
        status, output = commands.getstatusoutput(input)
        log.info(output)
        return status, output

    def cmd_ok(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertTrue(status == 0)

    def cmd_fails(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertFalse(status == 0)

# download ncbi taxonomy data and create a database if necessary; use
# this database for all non-destructive, non-modifying tests. For
# modifying tests, make a copy of the database.
ncbi_data, downloaded = ncbi.fetch_data(dest_dir=outputdir, clobber=False)
ncbi_master_db = path.join(outputdir, 'ncbi_master.db')
with ncbi.db_connect(ncbi_master_db, clobber = False) as con:
    log.info('using %s for all tests' % ncbi_master_db)
    ncbi.db_load(con, ncbi_data)

