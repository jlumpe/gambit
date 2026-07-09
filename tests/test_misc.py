"""Test miscellaneous stuff not tied to a specific module."""

from multiprocessing import cpu_count

import pytest
import numpy as np

from gambit._cython import threads


class TestCythonParallelism:
	"""Test that Cython modules are able to make use of parallelism.

	This can fail if the right compile and link args are not passed.
	"""

	def get_actual_nthreads(self, n: int = 100):
		"""Try to get the actual number of threads used in practice by OpenMP.

		Runs a multithreaded loop checks the number of unique thread IDs found.
		"""
		thread_ids = threads.get_thread_ids(n)
		ids_unique = set(thread_ids)
		nthreads = len(ids_unique)
		assert ids_unique == set(range(nthreads))
		return nthreads

	def test_default_nthreads(self):
		"""Test that the default maximum thread count is equal to the number of cores."""

		ncpus = cpu_count()
		assert threads.omp_get_max_threads() == ncpus
		assert self.get_actual_nthreads() == ncpus

	def test_omp_set_num_threads(self):
		"""Test the omp_set_num_threads() function."""

		nthreads_before = threads.omp_get_max_threads()
		nthreads = max(1, nthreads_before // 2)

		try:
			threads.omp_set_num_threads(nthreads)
			assert threads.omp_get_max_threads() == nthreads
			assert self.get_actual_nthreads() == nthreads

		finally:
			threads.omp_set_num_threads(nthreads_before)


def test_numpy_errors():
	"""Test the raise_numpy_errors fixture."""

	# Really just care about integer overflow
	with pytest.raises(FloatingPointError):
		np.uint8(255) + np.uint8(1)
	with pytest.raises(FloatingPointError):
		np.uint8(0) - np.uint8(1)


def test_version_attr():
	"""Test the gambit.__version__ attribute, kept for backwards compatibility.

	The version is now defined in pyproject.toml instead of a module-level variable, and is
	exposed at runtime through module-level __getattr__().
	"""

	from importlib.metadata import version, PackageNotFoundError
	from warnings import warn
	import gambit

	# Ensure we can still access standard attributes.
	assert gambit.__author__ == 'Jared Lumpe'

	# Ensure that the __getattr__() implementation still raises the correct error for missing attributes
	with pytest.raises(AttributeError):
		gambit.this_attribute_does_not_exist

	# Now try the version attribute
	try:
		expected_version = version('gambit')
	except PackageNotFoundError:
		warn('Unable to determine gambit package version.')
		assert gambit.__version__ is None
	else:
		assert gambit.__version__ == expected_version
