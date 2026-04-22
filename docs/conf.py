# -*- coding: utf-8 -*-
import sys
import os
import datetime

sys.path.insert(0, os.path.abspath('..'))

from taxtastic import __version__

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.todo', 'sphinx.ext.viewcode']

templates_path = ['_templates']
source_suffix = '.rst'
root_doc = 'index'

project = 'taxtastic'
author = 'Noah Hoffman, Erick Matsen, Brian Hodges, Connor McCoy, Chris Rosenthal'
copyright = '2011-{} {}'.format(datetime.date.today().strftime('%Y'), author)

version = __version__.split('+')[0]
release = version

exclude_patterns = ['_build', '.DS_Store']

html_theme = 'furo'
html_baseurl = 'https://fhcrc.github.io/taxtastic/'
pygments_style = 'friendly'
pygments_dark_style = 'monokai'

htmlhelp_basename = 'taxtasticdoc'

latex_documents = [
    ('index', 'taxtastic.tex', 'taxtastic Documentation',
     'Noah Hoffman, Erick Matsen, Brian Hodges', 'manual'),
]

man_pages = [
    ('index', 'taxtastic', 'taxtastic Documentation',
     ['Noah Hoffman, Erick Matsen, Brian Hodges, Frederick Ross'], 1)
]
