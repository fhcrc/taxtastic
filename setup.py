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

import versioneer

versioneer.versionfile_source = 'taxtastic/_version.py'
versioneer.versionfile_build = 'taxtastic/_version.py'
versioneer.tag_prefix = 'v' # tags are like v1.2.0
versioneer.parentdir_prefix = 'taxtastic-'

class run_audit(Command):
    """Audits source code using PyFlakes for following issues:
        - Names which are used but not defined or used before they are defined.
        - Names which are redefined without having been used.
    """
    description = "Audit source code with PyFlakes"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import os, sys
        try:
            import pyflakes.scripts.pyflakes as flakes
        except ImportError:
            print "Audit requires PyFlakes installed in your system."
            sys.exit(-1)

        warns = 0
        # Define top-level directories
        dirs = ['taxtastic']
        for dir in dirs:
            for root, _, files in os.walk(dir):
                for file in files:
                    if file != '__init__.py' and file.endswith('.py') :
                        warns += flakes.checkPath(os.path.join(root, file))
        if warns > 0:
            print "Audit finished with total %d warnings." % warns
        else:
            print "No problems found in sourcecode."

scripts = ['taxit']

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': find_packages(exclude=['tests']),
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': versioneer.get_version(),
          'cmdclass': versioneer.get_cmdclass(),
          'license': 'GPL',
          'classifiers': [
              'License :: OSI Approved :: GNU General Public License (GPL)',
              'Development Status :: 3 - Alpha',
              'Programming Language :: Python :: 2.7',
              'Topic :: Scientific/Engineering :: Bio-Informatics'],
          'download_url': 'https://github.com/fhcrc/taxtastic',
          'package_data': {'taxtastic': [join('data',f) for f in ['sha']]},
          'install_requires': [
              'sqlalchemy>=0.7',
              'decorator',
              'biopython',
              'xlrd']}

setup(**params)

