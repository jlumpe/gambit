"""setuptools installation script for gambit package"""

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


# Cython extensions
np_include = numpy.get_include()
extensions = [Extension(
	'gambit._cython.*',
	['src/gambit/_cython/*.pyx'],
	include_dirs=[np_include],
	extra_compile_args=['-fopenmp', '-Wno-sign-compare'],
	extra_link_args=['-fopenmp'],
)]
ext_modules = cythonize(
	extensions,
	compiler_directives=dict(
		language_level='3str',
		boundscheck=False,
		wraparound=False,
	),
)


setup(ext_modules=ext_modules)
