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
copyright = '2016 - 2026, Jared Lumpe'
author = 'Jared Lumpe'

# The full version, including alpha/beta/rc tags.
# Requires the package to be installed (Read the Docs does this, see .readthedocs.yml).
from importlib.metadata import version as _get_version
release = _get_version('gambit')
# Major.minor version
version = '.'.join(release.split('.')[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.autodoc',
	# 'sphinx.ext.doctest',
	'sphinx.ext.intersphinx',
	'sphinx.ext.todo',
	# 'sphinx.ext.coverage',
	# 'sphinx.ext.mathjax',
	# 'sphinx.ext.viewcode',
	'sphinx.ext.napoleon',
	'linuxdoc.rstFlatTable',  # For table cells spanning multiple rows
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# When debugging broken cross references using nitpick mode (-n option), ignore these errors.
# This mostly relates to external libraries that have not been linked to using intersphinx.
nitpick_ignore_regex = [
	('py:.*', r'click\..*'),
	('py:.*', r'sqlalchemy\..*'),
	('py:.*', r'h5py\..*'),
	('py:.*', r'scipy\..*'),
	# TypeVar / PEP 695 type parameters
	('py:.*', r'(.*\.)?T\d?'),
]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['css/custom.css']


# -- Extension options -------------------------------------------------------

autodoc_default_options = {
	'members': True,
	'show-inheritance': True,
}

autodoc_class_signature = 'separated'
autodoc_member_order = 'groupwise'
autodoc_typehints = 'description'

intersphinx_mapping = {
	'python': ('https://docs.python.org/3', None),
	'numpy': ('https://numpy.org/doc/stable/', None),
	'Bio': ('https://biopython.org/docs/latest/', None),
}

todo_include_todos = True


# -- Reference fixes ---------------------------------------------------------

# Targets whose runtime qualified name differs from where they are documented.
_REF_TARGET_FIXES = {
	'concurrent.futures._base.Executor': 'concurrent.futures.Executor',
}


def _fix_missing_reference(app, env, node, contnode):
	from sphinx.ext.intersphinx import missing_reference

	if node.get('refdomain') != 'py':
		return None

	target = node['reftarget']
	newnode = node.deepcopy()
	newnode['reftarget'] = _REF_TARGET_FIXES.get(target, target)

	# Numpy scalar types (e.g. numpy.float32) are documented as attributes, not classes
	if node['reftype'] == 'class' and target.startswith('numpy.'):
		newnode['reftype'] = 'attr'
	elif newnode['reftarget'] == target:
		return None

	return missing_reference(app, env, newnode, contnode)


def setup(app):
	app.connect('missing-reference', _fix_missing_reference)
