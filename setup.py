"""
Create unix package:    python setup.py sdist
Upload to pypi:         python setup.py sdist upload
"""

import subprocess
import os
from distutils.version import LooseVersion
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3

print('using {}, sqlite3 version {}'.format(sqlite3.__name__, sqlite3.sqlite_version))

min_sqlite3_version = '3.8.3'
if LooseVersion(sqlite3.sqlite_version) < LooseVersion(min_sqlite3_version):
    raise ImportError(('the sqlite3 library version for this python interpreter is '
                       '{}, but a version >= {} is required; '
                       'see https://github.com/fhcrc/taxtastic#installing').format(
                           sqlite3.sqlite_version, min_sqlite3_version))


subprocess.call(
    ('mkdir -p taxtastic/data && '
     'git describe --tags --dirty > taxtastic/data/ver.tmp'
     '&& mv taxtastic/data/ver.tmp taxtastic/data/ver '
     '|| rm -f taxtastic/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

# must import __version__ *after* running 'git describe' above
from taxtastic import __version__

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as fi:
    long_description = fi.read()


params = {'name': 'taxtastic',
          'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'maintainer': 'Chris Rosenthal',
          'maintainer_email': 'crosenth@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'long_description': long_description,
          'packages': find_packages(
              exclude=['tests', 'testfiles', 'devtools', 'docs']),
          'python_requires': '>=3.4',
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'license': 'GPL',
          'classifiers': [
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Development Status :: 3 - Alpha',
              'Programming Language :: Python :: 3.4',
              'Programming Language :: Python :: 3.5',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',
              'Topic :: Scientific/Engineering :: Bio-Informatics'],
          'download_url': 'https://github.com/fhcrc/taxtastic',
          'package_data': {
              'taxtastic': ['data/*']},
          'entry_points': {'console_scripts': ['taxit = taxtastic.scripts.taxit:main']},
          'test_suite': 'tests',
          'install_requires': [
              'DendroPy>=4.3.0',
              'PyYAML>=3.12',
              'decorator>=4.1.2',
              'fastalite>=0.3',
              'jinja2>=2.9',
              'psycopg2-binary>=2.7.3.1',
              'six',
              'sqlalchemy>=1.0',
          ]}

setup(**params)
