# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AIM Documentation'
copyright = '2025, Muhammad Hassan'
author = 'Muhammad Hassan'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',
              'sphinx.ext.mathjax']

templates_path = ['_templates']
exclude_patterns = []

mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
    }
}
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_theme = 'furo'
html_baseurl = "https://hassan-azizi.github.io/AIM/"
html_static_path = ['_static', 'images']
html_logo = 'images/AIM_logo.png'
# html_theme_options = { 'style_nav_header_background': '#BED2E0',
#                       'logo_only': True}
html_theme_options = {
    # 'style_nav_header_background': '#3CB371',
    'logo_only': True,                         
    'collapse_navigation': False,              
    'sticky_navigation': False,                 
    'navigation_depth': 3,                   
    'includehidden': True,                     
    'titles_only': True    
}