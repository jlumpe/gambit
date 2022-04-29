"""OpenMP stuff."""

from cython import parallel

import numpy as np
cimport numpy as np
cimport openmp


def omp_set_num_threads(n: int):
	"""Set maximum number of threads used by OpenMP.

	Just calls the ``omp_set_num_threads`` C function.
	"""
	if n <= 0:
		raise ValueError('Argument must be positive.')
	openmp.omp_set_num_threads(n)


def omp_get_max_threads():
	"""Get the maximum number of threads used by OpenMP.

	Just calls the ``omp_get_max_threads`` C function.
	"""
	return openmp.omp_get_max_threads()


def get_thread_ids(int num_threads):
	"""Run a multithreaded loop and get the thread ID running in each iteration."""

	cdef:
		np.ndarray[np.intp_t, ndim=1] thread_ids
		np.intp_t thread_id = -1
		int i

	thread_ids = np.full(num_threads, -1, dtype=np.intp)

	for i in parallel.prange(num_threads, nogil=True, schedule='static', chunksize=1):
		thread_id = parallel.threadid()
		thread_ids[i] = thread_id

	return thread_ids
