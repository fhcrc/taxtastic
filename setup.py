"""
Create unix package:    python setup.py sdist
Upload to pypi:         python setup.py sdist upload
"""

import subprocess
import os
from setuptools import setup, find_packages
# To use a consistent encoding
from os import path

subprocess.call(
    ('mkdir -p taxtastic/data && '
     'git describe --tags --dirty > taxtastic/data/ver.tmp'
     '&& mv taxtastic/data/ver.tmp taxtastic/data/ver '
     '|| rm -f taxtastic/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, 'w', encoding='utf-8'))

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
          'packages': find_packages(exclude=['testfiles', 'tests']),
          'package_dir': {'taxtastic': 'taxtastic'},
          'python_requires': '>=3.8',
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'license': 'GPL',
          'classifiers': [
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Development Status :: 3 - Alpha',
              'Programming Language :: Python :: 3.8',
              'Programming Language :: Python :: 3.9',
              'Programming Language :: Python :: 3.10',
              'Programming Language :: Python :: 3.11',
              'Programming Language :: Python :: 3.12',
              'Topic :: Scientific/Engineering :: Bio-Informatics'],
          'download_url': 'https://github.com/fhcrc/taxtastic',
          'include_package_data': True,
          'entry_points': {
              'console_scripts': ['taxit = taxtastic.scripts.taxit:main']},
          'test_suite': 'tests',
          'install_requires': [
              'decorator',
              'DendroPy',
              'fastalite',
              'jinja2',
              'psycopg-binary',
              'psycopg2-binary',
              'PyYAML',
              'sqlalchemy>=2',
              'sqlparse',
          ]}

setup(**params)
