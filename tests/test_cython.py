"""Tests to ensure that Cython modules have been compiled correctly."""

from multiprocessing import cpu_count

import pytest

from gambit._cython.test import get_thread_ids


def test_parallel():
	"""Test that Cython modules are able to make use of parallelism.

	This can fail if the right compile and link args are not passed.
	"""

	ncpus = cpu_count()

	# Run a multi-threaded loop and check the thread ID in each iteration
	thread_ids = get_thread_ids(ncpus)

	# Check that each loop iteration got its own thread
	assert set(thread_ids) == set(range(ncpus))
