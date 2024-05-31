# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import subprocess, os


subprocess.call('cd .. ; doxygen', shell=True)

project = 'CyRSoXS'
copyright = '2019-2022 Iowa State University'
author = 'Kumar Saurabh, Adarsh Krishnamurthy, Baskar Ganapathysubramanian, Eliot Gann, Dean M. Delongchamp, Peter J. Dudenas, Tyler B. Martin, Peter Beaucage'
release = '1.1.6.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

highlight_language = 'c++'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinxdoc'
html_static_path = ['_static']
html_extra_path = ['../html']
