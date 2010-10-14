import os
import logging

debuglevel = logging.WARNING

datadir = '../testfiles'
outputdir = '../test_output'
try:
    os.mkdir(outputdir)
except OSError:
    pass

