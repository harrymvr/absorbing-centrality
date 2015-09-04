# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import sys, os

sys.path.append(os.path.abspath('..'))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]
autosummary_generate = True

if os.getenv('SPELLCHECK'):
    extensions += 'sphinxcontrib.spelling',
    spelling_show_suggestions = True
    spelling_lang = 'en_US'

source_suffix = '.rst'
master_doc = 'index'
project = u'absorbing-centrality'
year = u'2015'
author = u'Charalampos Mavroforakis'
copyright = '2015 {0}'.format(author)
import pkg_resources
try:
    version = pkg_resources.get_distribution('absorbing_centrality').version
except pkg_resources.DistributionNotFound:
    print 'To build the documentation, the distribution information of'
    print 'absorbing_centrality has to be available.  Either install the'
    print 'package into your development environment or run "setup.py develop"'
    print 'to setup the metadata.  A virtualenv is recommended!'
    sys.exit(1)
del pkg_resources

# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

pygments_style = 'trac'
templates_path = ['.']
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = True
html_sidebars = {
   '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

intersphinx_mapping = {
    'python': ('http://docs.python.org/2.7', 'http://docs.python.org/objects.inv'),
    'networkx': ('http://networkx.github.io/documentation/latest','http://networkx.github.io/documentation/latest/objects.inv'),
    'numpy': ('http://docs.scipy.org/doc/numpy','http://docs.scipy.org/doc/numpy/objects.inv'),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference','http://docs.scipy.org/doc/scipy/reference/objects.inv'),
}

add_module_names = False
show_authors = False
