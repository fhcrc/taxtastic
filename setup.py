"""
Create unix package:    python setup.py sdist
"""

try:
    from setuptools import setup, find_packages
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages

import glob

from taxtastic import __version__

# all files with .py extension in top level are assumed to be scripts
scripts = ['taxit'] + list(set(glob.glob('*.py')) - set(['setup.py', 'distribute_setup.py']))

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': find_packages(exclude=['tests']),
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'install_requires': ['sqlalchemy', 'decorator']}

setup(**params)

