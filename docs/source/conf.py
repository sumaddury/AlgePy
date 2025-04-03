# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AlgePy'
copyright = '2025, Sucheer Maddury'
author = 'Sucheer Maddury'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
autodoc_member_order = 'bysource'

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

extensions = [
    'sphinx.ext.autodoc',  # Auto-generates docs from docstrings
    'sphinx.ext.napoleon',  # Supports Google-style docstrings
]