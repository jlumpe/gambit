"""Cython functions for calculating k-mer distance metrics"""

from cython.parallel import prange, parallel


def jaccard(COORDS_T[:] coords1, COORDS_T_2[:] coords2):
	"""Compute the Jaccard index between two k-mer sets in sparse coordinate format.

	Arguments are Numpy arrays containing k-mer indices in sorted order. Data types must be 16, 32,
	or 64-bit signed or unsigned integers, but do not need to match.

	This is by far the most efficient way to calculate the metric (this is a native function) and
	should be used wherever possible.

	Parameters
	----------
	coords1 : numpy.ndarray
		K-mer set in sparse coordinate format.
	coords2 : numpy.ndarray
		K-mer set in sparse coordinate format.

	Returns
	-------
	numpy.float32
		Jaccard index between the two sets, a real number between 0 and 1.

	See Also
	--------
	.jaccarddist
	"""
	return 1 - c_jaccarddist(coords1, coords2)


def jaccarddist(COORDS_T[:] coords1, COORDS_T_2[:] coords2):
	"""Compute the Jaccard distance between two k-mer sets in sparse coordinate format.

	The Jaccard distance is equal to one minus the Jaccard index.

	Arguments are Numpy arrays containing k-mer indices in sorted order. Data types must be 16, 32,
	or 64-bit signed or unsigned integers, but do not need to match.

	This is by far the most efficient way to calculate the metric (this is a native function) and
	should be used wherever possible.

	Parameters
	----------
	coords1 : numpy.ndarray
		K-mer set in sparse coordinate format.
	coords2 : numpy.ndarray
		K-mer set in sparse coordinate format.

	Returns
	-------
	numpy.float32
		Jaccard distance between the two sets, a real number between 0 and 1.

	See Also
	--------
	.jaccard
	"""
	return c_jaccarddist(coords1, coords2)


cdef SCORE_T c_jaccarddist(COORDS_T[:] coords1, COORDS_T_2[:] coords2) nogil:
	"""Compute the Jaccard distance between two k-mer sets in ordered coordinate format.

	Declared with nogil so it can be run in parallel.
	"""

	cdef:
		# Lengths of the two arrays
		intptr_t N = coords1.shape[0]
		intptr_t M = coords2.shape[0]

		# Index and value of items in each array as we are iterating
		intptr_t i = 0, j = 0
		COORDS_T a
		COORDS_T_2 b

		intptr_t u = 0  # Size of union

	# Iterate through both arrays simultaneously, advance index for the array
	# with the smaller value. Advance both if they are equal. Increment the
	# union count each loop.
	while i < N and j < M:
		a = coords1[i]
		b = coords2[j]

		u += 1

		if a <= b:
			i += 1

		if b <= a:
			j += 1

	# In most cases we won't have i == N and j == M at the end of the loop,
	# account for the items that we didn't get through
	u += N - i
	u += M - j

	# Avoid divide by zero, define distance between empty sets to be zero
	if u == 0:
		return 0

	# |A intersection B| = |A| + |B| - |A union B|
	# |A symdiff B|      = 2|A union B| - |A| - |B|
	return <SCORE_T>(2 * u - N - M) / u


def _jaccarddist_parallel(COORDS_T[:] query, COORDS_T_2[:] ref_coords, BOUNDS_T[:] ref_bounds, SCORE_T[:] out):
	"""Calculate Jaccard distances between a query k-mer set and a collection of reference sets.

	Data types of k-mer coordinate arrays may be 16, 32, or 64-bit signed or
	unsigned integers, but must match.

	Internally, releases the GIL in the main loop and calculates distances in parallel.

	Parameters
	----------
	query : numpy.ndarray
		Query k-mer set in sparse coordinate format.
	ref_coords : numpy.ndarray
		Reference k-mer sets in sparse coordinate format, concatenated into a single array.
	ref_bounds : numpy.ndarray
		Bounds of individual k-mer sets within the ``ref_coords`` array. The ``n``\ th k-mer set is
		the slice of ``ref_coords`` between ``ref_bounds[n]`` and ``ref_bounds[n + 1]``. Length must
		be one greater than that of``ref_coords``.
	out : numpy.ndarray
		Pre-allocated array to write distances to.
	"""
	cdef intptr_t N = ref_bounds.shape[0] - 1
	cdef BOUNDS_T begin, end
	cdef int i

	for i in prange(N, nogil=True, schedule='dynamic'):
		begin = ref_bounds[i]
		end = ref_bounds[i+1]
		out[i] = c_jaccarddist(query, ref_coords[begin:end])
