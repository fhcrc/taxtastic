"""
Create unix package:    python setup.py sdist
"""

from distutils.core import setup
import glob

from taxtastic import __version__

# all files with .py extension in top level are assumed to be scripts
scripts = ['taxit'] + list(set(glob.glob('*.py')) - set(['setup.py']))

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': ['taxtastic', 'taxtastic.scripts', 'taxtastic.subcommands'],
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'requires': ['Python (>= 2.7)', 'sqlalchemy', 'decorator']}

setup(**params)

