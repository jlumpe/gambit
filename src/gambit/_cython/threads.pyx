"""OpenMP stuff."""

from cython import parallel
import array

cimport cython
from cpython cimport array
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


@cython.boundscheck(True)
def get_thread_ids(int n):
	"""Run a multithreaded loop and get the thread ID running in each iteration.

	Used to check that Cython code parallelization is working correctly. Result should contain
	integers from 0 to ``num_threads``, repeated up to length ``n``.

	Parameters
	----------
	n: int
		Size of loop. Make this at least as large as the expected number of threads.

	Returns
	-------
	array.array
		Array of size ``n`` containing the thread ID running in each loop iteration.
	"""

	cdef:
		array.array thread_ids_arr = array.array('i')
		int[:] thread_ids
		int i

	for i in range(n):
		thread_ids_arr.append(-1)

	thread_ids = thread_ids_arr

	for i in parallel.prange(n, nogil=True, schedule='static', chunksize=1):
		thread_ids[i] = parallel.threadid()

	return thread_ids
