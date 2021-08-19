"""Test miscellaneous stuff not tied to a specific module."""

from multiprocessing import cpu_count

import pytest
import numpy as np

from gambit._cython.test import get_thread_ids


def test_cython_parallel():
	"""Test that Cython modules are able to make use of parallelism.

	This can fail if the right compile and link args are not passed.
	"""

	ncpus = cpu_count()

	# Run a multi-threaded loop and check the thread ID in each iteration
	thread_ids = get_thread_ids(ncpus)

	# Check that each loop iteration got its own thread
	assert set(thread_ids) == set(range(ncpus))


def test_numpy_errors():
	"""Test the raise_numpy_errors fixture."""

	# Really just care about integer overflow
	with pytest.raises(FloatingPointError):
		np.uint8(255) + np.uint8(1)
	with pytest.raises(FloatingPointError):
		np.uint8(0) - np.uint8(1)
