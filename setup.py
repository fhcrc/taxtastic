"""
Create unix package:    python setup.py sdist
Upload to pypi:         python setup.py sdist upload
"""

from os.path import join
import glob

try:
    from setuptools import setup, find_packages, Command
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages, Command


def get_git_version():
    import subprocess
    git_version = subprocess.check_output(['git', 'describe']).rstrip()
    return git_version

class git_version(Command):
    description = "Show package version from git"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        git_version = get_git_version()
        print 'Version:', git_version

scripts = ['taxit']

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': find_packages(exclude=['tests']),
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': '0.4.0',
          'package_data': {'taxtastic': [join('data',f) for f in ['sha']]},
          'cmdclass': {'git_version': git_version},
          'install_requires': ['sqlalchemy', 'decorator', 'biopython']}

setup(**params)

