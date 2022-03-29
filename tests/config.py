import contextlib
from io import StringIO
import os
from os import path
import sys
import logging
import re
import unittest
import shutil
import tempfile
import subprocess


from taxtastic.scripts.taxit import main


log = logging


def funcname(idstr):
    return '.'.join(idstr.split('.')[1:])


# set verbosity of logging output
try:
    logflag = re.findall(r'-[vq]+\b', ' '.join(sys.argv[1:]))[0]
except IndexError:
    logflag = ''

logging.basicConfig(
    stream=sys.stdout,
    format='%(levelname)s %(module)s %(lineno)s %(message)s'
    if logflag.startswith('-v') else '%(message)s',
    level={'': logging.CRITICAL,
           '-q': logging.CRITICAL,
           '-v': logging.WARNING,
           '-vv': logging.INFO,
           '-vvv': logging.DEBUG}.get(logflag)
)

# module data
datadir = path.abspath(path.join(path.dirname(__file__), '..', 'testfiles'))
outputdir = path.abspath(
    path.join(path.dirname(__file__), '..', 'test_output'))


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
    except OSError as msg:
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
        self.silent = '-v' not in sys.argv

        if self.silent:
            self.old_stdout = sys.stdout
            self.old_stderr = sys.stderr
            sys.stdout = StringIO()
            sys.stderr = StringIO()

    def tearDown(self):
        if self.silent:
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

    def suppress_stdout(self):
        self.old_stdout = sys.stdout
        sys.stdout = StringIO()

    def suppress_stderr(self):
        self.old_stderr = sys.stderr
        sys.stderr = StringIO()

    def tearDown(self):
        if hasattr(self, 'old_stdout'):
            sys.stdout = self.old_stdout
        if hasattr(self, 'old_stderr'):
            sys.stdout = self.old_stderr


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
    openargs = {'mode': 'r'} if sys.version_info.major == 2 else {'newline': None}

    def __getitem__(self, i):
        """
        Enables string formatting, eg:
        print 'current function: %(funcname)s' % self
        """
        return getattr(self, i)

    def wrap_cmd(self, cmd):
        cmd = cmd % self
        try:
            retval = main(cmd.split())
        except SystemExit as err:
            return None, err
        else:
            return retval, None

    def cmd_ok(self, cmd):
        retval, err = self.wrap_cmd(cmd)
        if err:
            self.assertEqual(errval(err), 0)
        else:
            self.assertFalse(bool(retval))

    def cmd_fails(self, cmd):
        retval, err = self.wrap_cmd(cmd)
        if err:
            self.assertNotEqual(errval(err), 0)
        else:
            self.assertTrue(bool(retval))


def errval(err):
    if sys.version_info.major == 2:
        return err.message
    else:
        return err.args[0]


# Small NCBI taxonomy database
# this database for all non-destructive, non-modifying tests. For
# modifying tests, make a copy of the database.
ncbi_master_db = data_path('small_taxonomy.db')
ncbi_data = data_path('taxdmp.zip')
