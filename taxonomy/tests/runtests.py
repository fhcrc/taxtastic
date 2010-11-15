#!/Users/nhoffman/devel/bin/python

import cProfile
import glob
import logging
import os
import pstats
import re
import sys
import unittest
import pprint
from optparse import OptionParser

sys.path.insert(0,'.')

log = logging
    
def main():
    """
    Run all unit tests for a python module.
    """

    usage = "%prog [options] [data files]"
    usage += main.__doc__.strip()
    
    parser = OptionParser(usage=main.__doc__, version="$Id$")
    parser.set_defaults(
        verbose=0,
        pattern='*_test.py',
        profile=False)

    parser.add_option("-p", "--pattern",
        action="store", dest="pattern", 
        help="[*_test.py] pattern used to identify unit test modules.")
    parser.add_option("-t", "--target",
        action="store", dest="target", metavar='NAME',
        help="Run tests contained in NAME, provided in dotted format (eg, ClassName or ClassName.methodName).")
    parser.add_option("-P", "--profile",
        action="store_true", dest="profile", 
        help="Run profiler on indicated unit tests")        
    parser.add_option("-v", "--verbose",
        action="count", dest="verbose", 
        help="increase verbosity of screen output")
        
    (options, args) = parser.parse_args()
    
    loglevel = {0:logging.WARNING, 
        1:logging.INFO, 
        2:logging.DEBUG}.get(options.verbose, logging.DEBUG)
    
    verbose_format = '%(levelname)s %(module)s %(lineno)s %(message)s'
    logformat = {0:'%(message)s', 
        1:verbose_format, 
        2:verbose_format}.get(options.verbose, verbose_format)
    
    # set up logging
    logging.basicConfig(file=sys.stdout,
        format=logformat,
        level=loglevel)
    
    target = options.target
    
    if args:
        fnames = args
    else:
        fnames = glob.glob(options.pattern)
    
    log.warning('input files: %s' % ', '.join(fnames))
    
    suites = []
    for fname in fnames:
        module = os.path.splitext(fname)[0] 
        if target:
            load_name = '%(module)s.%(target)s' % locals()
        else:
            load_name = module
        
        log.info('loading %s' % load_name)
            
        try:
            suites.append(unittest.TestLoader().loadTestsFromName(load_name))
        except AttributeError, msg:
            log.error(msg)
            continue 
    
    log.info('test suites: %s' % pprint.pformat(suites))
    
    if options.profile:
        profname = 'profile'
        cProfile.runctx('unittest.TextTestRunner(\
        verbosity=2).run(unittest.TestSuite(suites))', 
        globals(), locals(), profname)
        p = pstats.Stats(profname).strip_dirs()
        p.sort_stats('cumulative').print_stats(15)    
    else:
        unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suites))   
    
if __name__ == '__main__':
    main()