import os
from os import path
import sys
import logging
import re
import unittest
import commands
import shutil

log = logging

def funcname(idstr):
    return '.'.join(idstr.split('.')[1:])

def mkdir(dirpath):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created.
    """
    
    try:
        os.mkdir(dirpath)
    except OSError, msg:
        log.warning(msg)

    if not path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath
    
def rmdir(dirpath):
    """
    Remove a (potentially missing) directory without errors. Raise
    OSError if directory can't be removed.
    """

    try:
        shutil.rmtree(dirpath)
    except OSError, msg:
        log.warning(msg)

    if path.exists(dirpath):
        raise OSError('Failed to remove %s' % dirpath)
    
class TestScriptBase(unittest.TestCase):

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
    outputdir = None
    
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

    def setUp(self):
        self.funcname = '.'.join(self.id().split('.')[1:])
        self.package = path.join(self.outputdir, self.funcname)

    def tearDown(self):
        pass


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
datadir = 'testfiles'
outputdir = 'test_output'

mkdir(outputdir)

