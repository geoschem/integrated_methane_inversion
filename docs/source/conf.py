# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Custom CSS for tables ---------------------------------------------------

def setup(app):
   app.add_css_file('css/custom.css')

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

master_doc = 'index'

# -- Project information -----------------------------------------------------

project = 'IMI'
copyright = '2021, TODO'
author = 'TODO'

# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx_rtd_theme",   
    "sphinxcontrib.bibtex",
    "recommonmark",
]
bibtex_default_style = 'gcrefstyle'

from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.names.lastfirst import NameStyle as LastFirst
from pybtex.style.template import join, words, optional, sentence
from pybtex.style.labels import BaseLabelStyle

class GCLabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key.replace("_", " ").replace("et al.", "et al.,")

class GCRefStyle(UnsrtStyle):
    default_name_style = LastFirst
    default_sort_style = None
    default_label_style = GCLabelStyle

    def __init__(self):
       super().__init__()
       self.abbreviate_names = True
      #  self.label_style = KeyLabelStyle()
      #  self.format_labels = self.label_style.format_labels

    def format_web_refs(self, e):
       return sentence[ optional[ self.format_doi(e) ], ]

from pybtex.plugin import register_plugin
register_plugin('pybtex.style.formatting', 'gcrefstyle', GCRefStyle)


bibtex_bibliography_header = ".. rubric:: References"
bibtex_footbibliography_header = bibtex_bibliography_header

bibtex_bibfiles = ['geos-chem-shared-docs/geos-chem.bib']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
#html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add favicon (small image in browser tabs)
html_favicon = '_static/favicon.ico'

# Add IMI logo at top left
html_logo = "_static/imi_logo_small.png"

# Theme customization: show logo on white background, no version number
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'style_nav_header_background': '#FFFFFF',
}
