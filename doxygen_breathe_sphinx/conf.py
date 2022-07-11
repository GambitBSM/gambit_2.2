# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'GAMBIT'
copyright = '2022, The GAMBIT Community'
author = 'The GAMBIT Community'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ "breathe" ]

breathe_projects = {"GAMBIT": "./xml"}

breathe_default_project = "GAMBIT"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add custom template
templates_path = ['_templates']

# Change name of main html file
root_doc = 'index'

# Try to clean up output
html_show_copyright = False
html_show_sphinx = False
html_permalinks = True
html_sidebars = {'**': ['searchbox.html', 'globaltoc.html', 'localtoc.html']}
html_copy_source = False

# RTD theme options
html_theme_options = {
    'prev_next_buttons_location': 'none',
}