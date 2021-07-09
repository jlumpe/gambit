"""setuptools installation script for gambit package"""

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


# Cython extensions
np_include = numpy.get_include()
extensions = [Extension(
	'gambit._cython.*',
	['gambit/_cython/*.pyx'],
	include_dirs=[np_include],
	extra_compile_args=['-fopenmp', '-Wno-sign-compare'],
	extra_link_args=['-fopenmp'],
)]


setup(ext_modules=cythonize(extensions))
