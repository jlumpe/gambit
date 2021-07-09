"""Tests for Cython-compiled code."""

from cython import parallel
from libc.stdio cimport printf

import numpy as np
cimport numpy as np


def print_thread_ids(int num_threads=10):
	"""Start a number of threads and have each print its ID."""
	cdef int thread_id = -1
	with nogil, parallel.parallel(num_threads=num_threads):
		thread_id = parallel.threadid()
		printf("Thread ID: %d\n", thread_id)


def get_thread_ids(int num_threads=10):
	"""Run a multithreaded loop and get the thread ID that read each loop iteration."""

	cdef:
		np.ndarray[np.intp_t, ndim=1] thread_ids
		np.intp_t thread_id = -1
		int i

	thread_ids = np.full(num_threads, -1, dtype=np.intp)

	for i in parallel.prange(num_threads, nogil=True, schedule='static', chunksize=1):
		thread_id = parallel.threadid()
		thread_ids[i] = thread_id

	return thread_ids
