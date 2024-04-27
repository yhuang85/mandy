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

import string

# -- Project information -----------------------------------------------------

project = 'Mandy'
copyright = '2021, Mandy'
author = 'Mandy'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.mathjax',
	'sphinx.ext.todo',
	'sphinx_design',
]

todo_include_todos = True

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]



html_logo = '_static/dms.png'
html_title = "Mandy's wonderland"
html_theme_options = {
	"icon_links": [
		{
			"name": "GitHub",
			"url": "https://github.com/yhuang85/mandy",
			"icon": "fab fa-github-square",
		}
	],
	"navigation_depth": 4,
	"show_prev_next": False,
	"use_edit_page_button": True,
}

html_context = {
	"github_user": "yhuang85",
	"github_repo": "mandy",
	"github_version": "main",
	"doc_path": "build-source/source",
}

html_sidebars = {
	# Hide primary sidebar
	"**": []
}

mathjax3_config = {
	'tex': {
		'macros': {
			'op': ['\\operatorname {#1}', 1],
			'p': '{\\partial}',
			'ps': ['\\prescript {#1}{#2}{#3}', 3],
			'psup': ['\\prescript {#1}{}{#2}', 2],
			'blue': ['{\\color{blue} {#1}}', 1],
			**{w: f'{{\\operatorname{{{w}}}}}' for w in ['dist', 'ind', 'std']},
			**{f'{w}bb': f'{{\\mathbb {w}}}' for w in string.ascii_letters},
			**{f'{w}bf': f'{{\\mathbf {w}}}' for w in string.ascii_letters},
			**{f'{w}cal': f'{{\\mathcal {w}}}' for w in string.ascii_letters},
			**{f'{w}frak': f'{{\\mathfrak {w}}}' for w in string.ascii_letters},
            **{f'{w}scr': f'{{\\mathscr {w}}}' for w in string.ascii_letters},
		},
		'packages': {
			'[+]': [
				'ams',
				'color',
				'mathtools',
				'unicode',
			],
		},
		'tags': 'ams',
	},
	'loader': {
		'load': [
			'[tex]/ams',
			'[tex]/color',
			'[tex]/mathtools',
			'[tex]/unicode',
		],
	}
}
