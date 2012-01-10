import contextlib
from cStringIO import StringIO
import os
from os import path, rmdir
import sys
import logging
import re
import unittest
import commands
import shutil
import tempfile

log = logging

def funcname(idstr):
    return '.'.join(idstr.split('.')[1:])

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
datadir = path.abspath(path.join(path.dirname(__file__), '..', 'testfiles'))
outputdir = path.abspath(path.join(path.dirname(__file__), '..', 'test_output'))

def mkdir(dirpath, clobber=False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)

    try:
        os.mkdir(dirpath)
    except OSError, msg:
        log.debug(msg)

    if not path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath

if not os.path.isdir(outputdir):
    mkdir(outputdir)

def data_path(*args):
    return os.path.join(datadir, *args)

def output_path(*args):
    return os.path.join(outputdir, *args)

@contextlib.contextmanager
def tempdir(*args, **kwargs):
    try:
        d = tempfile.mkdtemp(*args, **kwargs)
        yield d
    finally:
        shutil.rmtree(d)

class OutputRedirectMixin(object):
    def setUp(self):
        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr
        sys.stdout = StringIO()
        sys.stderr = StringIO()

    def tearDown(self):
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr


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

    def mkoutdir(self, clobber=True):
        """
        Create output directory (destructively if clobber is True)
        """

        outdir = self.outdir()
        mkdir(outdir, clobber)
        return outdir


class TestScriptBase(OutputRedirectMixin, TestBase):

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
        log.info('--> '+ input)
        status, output = commands.getstatusoutput(input)
        log.info(output)
        return status, output

    def cmd_ok(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertTrue(status == 0)

    def cmd_fails(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertFalse(status == 0)

# Small NCBI taxonomy database
# this database for all non-destructive, non-modifying tests. For
# modifying tests, make a copy of the database.
ncbi_master_db = data_path('small_taxonomy.db')
ncbi_data = data_path('taxdmp.zip')
