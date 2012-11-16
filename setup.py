"""
Create unix package:    python setup.py sdist
Upload to pypi:         python setup.py sdist upload
"""

from os.path import join

try:
    from setuptools import setup, find_packages, Command
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages, Command

import taxtastic

scripts = ['taxit']

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': find_packages(exclude=['tests']),
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': taxtastic.__version__,
          'license': 'GPL',
          'classifiers': [
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Development Status :: 3 - Alpha',
              'Programming Language :: Python :: 2.7',
              'Topic :: Scientific/Engineering :: Bio-Informatics'],
          'download_url': 'https://github.com/fhcrc/taxtastic',
          'package_data': {'taxtastic': [join('data',f) for f in ['sha']]},
          'install_requires': ['sqlalchemy>=0.7', 'decorator', 'biopython']}

setup(**params)

