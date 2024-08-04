"""Calculate the Jaccard index/distance between sets."""

from collections.abc import Set
from typing import Iterable, Sequence, Optional

import numpy as np

import gambit._cython.metric as _cmetric
from gambit.sigs.base import KmerSignature, SignatureArray, AbstractSignatureArray, SignatureList, \
	BOUNDS_DTYPE
from gambit.util.misc import chunk_slices
from gambit.util.progress import get_progress


#: Numpy dtype for output of Cython Jaccard distance calculation code
# Equivalent to SCORE_T in types.pxd
SCORE_DTYPE = np.dtype(np.float32)


_COORDS_UNSIGNED_DTYPES = [np.dtype(f'u{s}') for s in [2, 4, 8]]
_COORDS_SIGNED_DTYPES = [np.dtype(f'i{s}') for s in [2, 4, 8]]


def _cast_sigs_array(arr: np.ndarray) -> np.ndarray:
	"""Convert signature array to proper data type for Cython metric code.

	Cython code accepts k-mer coordinate arrays in 16, 32, or 64-bit unsigned data types, these are
	returned as-is. Equivalent signed data types can safely be casted (as the values should all be
	non-negative), for these a view into the array with unsigned data type is returned (no coyping).
	All other data types result in a ValueError.
	"""

	dt = arr.dtype
	if dt in _COORDS_UNSIGNED_DTYPES:
		return arr
	if dt in _COORDS_SIGNED_DTYPES:
		new_dt = np.dtype(f'u{dt.itemsize}')
		return arr.view(new_dt)
	raise ValueError(f'Invalid dtype for k-mer coordinate array: {dt.str}')


def jaccard(coords1: np.ndarray, coords2: np.ndarray) -> np.float32:
	"""Compute the Jaccard index between two k-mer sets in sparse coordinate format.

	Arguments are Numpy arrays containing k-mer indices in sorted order. Data types must be 16, 32,
	or 64-bit signed or unsigned integers, but do not need to match.

	This is by far the most efficient way to calculate the metric (this is a native function) and
	should be used wherever possible.

	Parameters
	----------
	coords1
		K-mer set in sparse coordinate format.
	coords2
		K-mer set in sparse coordinate format.

	Returns
	-------
	numpy.float32
		Jaccard index between the two sets, a real number between 0 and 1.

	See Also
	--------
	.jaccarddist
	"""
	coords1 = _cast_sigs_array(coords1)
	coords2 = _cast_sigs_array(coords2)
	return _cmetric.jaccard(coords1, coords2)


def jaccarddist(coords1: np.ndarray, coords2: np.ndarray):
	"""Compute the Jaccard distance between two k-mer sets in sparse coordinate format.

	The Jaccard distance is equal to one minus the Jaccard index.

	Arguments are Numpy arrays containing k-mer indices in sorted order. Data types must be 16, 32,
	or 64-bit signed or unsigned integers, but do not need to match.

	This is by far the most efficient way to calculate the metric (this is a native function) and
	should be used wherever possible.

	Parameters
	----------
	coords1
		K-mer set in sparse coordinate format.
	coords2
		K-mer set in sparse coordinate format.

	Returns
	-------
	numpy.float32
		Jaccard distance between the two sets, a real number between 0 and 1.

	See Also
	--------
	.jaccard
	"""
	coords1 = _cast_sigs_array(coords1)
	coords2 = _cast_sigs_array(coords2)
	return _cmetric.jaccarddist(coords1, coords2)



def jaccard_generic(set1: Iterable, set2: Iterable) -> float:
	"""Get the Jaccard index of of two arbitrary sets.

	This is primarily used as a slow, pure-Python alternative to :func:`.jaccard` to be used
	for testing, but can also be used as a generic way to calculate the Jaccard index which works
	with any collection or element type.

	See Also
	--------
	.jaccard
	.jaccard_bits
	"""
	if not isinstance(set1, Set):
		set1 = set(set1)

	n1 = len(set1)
	n2 = len(set2)
	intersection = len(set1.intersection(set2))
	union = n1 + n2 - intersection

	return 1. if union == 0 else intersection / union


def jaccard_bits(bits1: np.ndarray, bits2: np.ndarray) -> float:
	"""Calculate the Jaccard index between two sets represented as bit arrays ("dense" format for k-mer sets).

	See Also
	--------
	.jaccard
	"""
	n1 = np.count_nonzero(bits1)
	n2 = np.count_nonzero(bits2)
	intersection = np.count_nonzero(bits1 & bits2)
	union = n1 + n2 - intersection
	return 1. if union == 0 else intersection / union


def jaccarddist_array(query: KmerSignature, refs: Sequence[KmerSignature], out: Optional[np.ndarray] = None) -> np.ndarray:
	"""
	Calculate Jaccard distances between a query k-mer signature and a list of reference signatures.

	For enhanced performance ``refs`` should be an instance of
	:class:`gambit.sigs.base.SignatureArray`. This allows use of optimized Cython code that runs
	in parallel over all signatures in ``refs``. In that case, because of Cython limitations
	``refs.bounds.dtype`` must be ``np.intp``, which is usually a 64-bit signed integer. If it is
	not it will be converted automatically.

	Parameters
	----------
	query
		Query k-mer signature in sparse coordinate format (sorted array of k-mer indices).
	refs
		List of reference signatures.
	out
		Optional pre-allocated array to write results to. Should be the same length as ``refs``
		with dtype ``np.float32``.

	Returns
	-------
	numpy.ndarray
		Jaccard distance for ``query`` against each element of ``refs``.

	See Also
	--------
	.jaccarddist
	.jaccarddist_matrix
	"""
	query = _cast_sigs_array(query)

	if out is None:
		out = np.empty(len(refs), SCORE_DTYPE)
	elif out.shape != (len(refs),):
		raise ValueError('Output array length must match signature array.')
	elif out.dtype != SCORE_DTYPE:
		raise ValueError(f'Output array dtype must be {SCORE_DTYPE}, got {out.dtype}')

	if isinstance(refs, SignatureArray):
		values = _cast_sigs_array(refs.values)
		bounds = refs.bounds.astype(BOUNDS_DTYPE, copy=False)

		_cmetric._jaccarddist_parallel(query, values, bounds, out)

	else:
		for i, ref in enumerate(refs):
			ref = _cast_sigs_array(ref)
			out[i] = _cmetric.jaccarddist(query, ref)

	return out


def jaccarddist_matrix(queries: Sequence[KmerSignature],
                       refs: Sequence[KmerSignature],
                       ref_indices: Optional[Sequence[int]] = None,
                       out: Optional[np.ndarray] = None,
                       chunksize: Optional[int] = None,
                       progress = None,
                       ) -> np.ndarray:
	"""
	Calculate a Jaccard distance matrix between a list of query signatures and a list of
	reference signatures.

	This function improves querying performance when the reference signatures are stored in a file
	(e.g. using :class:`gambit.sigs.hdf5.HDF5Signatures`) by loading them in chunks (via the
	``chunksize`` parameter) instead of all in one go.

	Performance is greatly improved if ``refs`` is a type that yields instances of
	``SignatureArray`` when indexed with a slice object (``SignatureArray`` or
	``HDF5Signatures``), see :func:`.jaccarddist_array`. There is no such dependence on the type of
	``queries``, which can be a simple list.

	Parameters
	----------
	queries
		Query signatures in sparse coordinate format.
	refs
		Reference signatures in sparse coordinate format.
	ref_indices
		Optional, indices of ``refs`` to use.
	out
		(Optional) pre-allocated array to write output to.
	chunksize
		Divide ``refs`` into chunks of this size.
	progress
		Display a progress meter of the number of elements of the output array calculated so far.
		See :func:`gambit.util.progress.get_progress` for a description of allowed values.

	Returns
	-------
	numpy.ndarray
		Matrix of distances between query signatures in rows and reference signatures in columns.

	See Also
	--------
	.jaccarddist
	.jaccarddist_array
	.jaccarddist_pairwise
	"""
	nqueries = len(queries)
	nrefs = len(refs) if ref_indices is None else len(ref_indices)

	if not isinstance(refs, AbstractSignatureArray):
		refs = SignatureList(refs)  # To support advanced indexing

	if out is None:
		out = np.empty((nqueries, nrefs), SCORE_DTYPE)
	elif out.shape != (nqueries, nrefs):
		raise ValueError('Output array must have shape (nqueries, nrefs).')
	elif out.dtype != SCORE_DTYPE:
		raise ValueError(f'Output array dtype must be {SCORE_DTYPE}, got {out.dtype}')

	if chunksize is None:
		ref_slices = [slice(0, nrefs)]
	else:
		ref_slices = list(chunk_slices(nrefs, chunksize))

	with get_progress(progress, nqueries * nrefs) as meter:
		for ref_slice in ref_slices:
			idx = ref_slice if ref_indices is None else ref_indices[ref_slice]
			ref_chunk = refs[idx]

			for (i, query) in enumerate(queries):
				jaccarddist_array(query, ref_chunk, out=out[i, ref_slice])
				meter.increment(len(ref_chunk))

	return out


def num_pairs(n: int) -> int:
	"""Get the number of distinct (unordered) pairs of ``n`` objects."""
	return n * (n - 1) // 2


def jaccarddist_pairwise(sigs: Sequence[KmerSignature],
                         indices: Optional[Sequence[int]] = None,
                         flat: bool = False,
                         out: Optional[np.ndarray] = None,
                         progress = None,
                         ) -> np.ndarray:
	"""
	Calculate all pairwise Jaccard distances for a list of signatures.

	This should be roughly twice as fast as calling :func:`.jaccarddist_matrix` with the same array
	for the first and second arguments, because each pairwise distance is computed once instead of
	twice.

	For optimal performance the type of ``sigs`` is subject to the same requirements as
	:func:`.jaccarddist_array` and :func:`.jaccarddist_matrix`.

	Parameters
	----------
	sigs
		List of signatures in sparse coordinate format.
	indices
		Optional, indices of ``sigs`` to use.
	flat
		If True the output is a non-redundant flat (1D) array with exactly one element per pair of
		signatures. This format can be converted to/from the equivalent full distance matrix with
		:func:`scipy.spatial.distance.squareform`.
	out
		(Optional) pre-allocated array to write output to.
	progress
		Display a progress meter of the number of elements of the output array calculated so far.
		See :func:`gambit.util.progress.get_progress` for a description of allowed values.

	Returns
	-------
	numpy.ndarray
		Pairwise distances in matrix (if ``flat=False``) or condensed (``flat=True``) format.

	See Also
	--------
	.jaccarddist_matrix
	"""
	if not isinstance(sigs, AbstractSignatureArray):
		sigs = SignatureList(sigs)  # To support advanced indexing

	if indices is not None:
		indices = np.asarray(indices)

	n = len(sigs) if indices is None else len(indices)
	npairs = num_pairs(n)

	out_shape = (npairs,) if flat else (n, n)
	if out is None:
		out = np.empty(out_shape, SCORE_DTYPE)

	else:
		if out.shape != out_shape:
			raise ValueError(f'Expected output array to have size {out_shape} for {n} signatures and flat={flat}')

		if out.dtype != SCORE_DTYPE:
			raise ValueError(f'Output array dtype must be {SCORE_DTYPE}, got {out.dtype}')

	if flat:
		next_out = 0
	else:
		np.fill_diagonal(out, 0)

	with get_progress(progress, npairs) as meter:
		for i in range(n - 1):
			row_sig = sigs[i] if indices is None else sigs[indices[i]]

			cols = slice(i + 1, n)
			ncol = n - i - 1
			col_sigs = sigs[cols] if indices is None else sigs[indices[cols]]

			row_out = out[next_out:next_out+ncol] if flat else out[i, cols]

			jaccarddist_array(row_sig, col_sigs, out=row_out)
			meter.increment(ncol)

			if flat:
				next_out += ncol
			else:
				# Copy to other side of diagonal
				out[cols, i] = out[i, cols]

	return out
