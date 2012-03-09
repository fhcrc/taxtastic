"""
Create unix package:    python setup.py sdist
"""

import os
import subprocess
import shutil
from os.path import join
import glob

try:
    from setuptools import setup, find_packages
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages

# The code contained by the with statement below configures and
# updates the insertion of the abbrevated git sha hash in the package
# version number. This assumes that the repository root contains a
# '.gitattributes' file with the contents
# 
# taxtastic/data/sha    filter=sha
#
# and that the (empty) file 'taxtastic/data/sha' exists.
subprocess.call('git log --pretty=format:%h -n 1 > taxtastic/data/sha', shell=True)
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
          'package_data': {'taxtastic': [join('data',f) for f in ['sha']]},
          'install_requires': ['sqlalchemy', 'decorator']}

setup(**params)

