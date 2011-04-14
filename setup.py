"""
Create unix package:    python setup.py sdist
"""

from distutils.core import setup
import os
import sys
import glob

from taxtastic.__init__ import __version__

# all files with .py extension in top level are assumed to be scripts
scripts = ['taxit'] + list(set(glob.glob('*.py')) - set(['setup.py']))
               
params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'package_dir': {'taxtastic': 'taxtastic'},
          'packages': ['taxtastic'],
          'scripts': scripts,
          # 'package_data':{'taxonomy': glob.glob('data/*')},
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'requires': ['Python (>= 2.7)', 'sqlalchemy']}

setup(**params)

