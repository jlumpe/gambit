"""Calculate the Jaccard index/distance between sets."""

from collections.abc import Set
from typing import Iterable, Sequence, Optional

import numpy as np

from gambit._cython.metric import BOUNDS_DTYPE, SCORE_DTYPE, jaccard, jaccarddist, \
	_jaccarddist_parallel
from gambit.sigs.base import KmerSignature, SignatureArray, AbstractSignatureArray, SignatureList
from gambit.util.misc import chunk_slices
from gambit.util.progress import get_progress


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


def jaccarddist_array(query: KmerSignature, refs: Sequence[KmerSignature], out: np.ndarray = None) -> np.ndarray:
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
	if out is None:
		out = np.empty(len(refs), SCORE_DTYPE)
	elif out.shape != (len(refs),):
		raise ValueError('Output array length must match signature array.')
	elif out.dtype != SCORE_DTYPE:
		raise ValueError(f'Output array dtype must be {SCORE_DTYPE}, got {out.dtype}')

	if isinstance(refs, SignatureArray):
		values = refs.values
		bounds = refs.bounds.astype(BOUNDS_DTYPE, copy=False)

		_jaccarddist_parallel(query, values, bounds, out)

	else:
		for i, ref in enumerate(refs):
			out[i] = jaccarddist(query, ref)

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
	``HDF5Signatures``), see :meth:`.jaccarddist_array`. There is no such dependence on the type of
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
	np.ndarray
		Matrix of distances between query signatures in rows and reference signatures in columns.

	See Also
	--------
	.jaccarddist
	.jaccarddist_array
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
