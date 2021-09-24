"""Calculate the Jaccard index/distance between sets."""

from collections.abc import Set
from typing import Iterable, Sequence, Optional

import numpy as np

from gambit._cython.metric import BOUNDS_DTYPE, SCORE_DTYPE, jaccard_sparse, jaccarddist_sparse, \
	_jaccard_sparse_parallel
from gambit.signatures import KmerSignature, SignatureArray
from gambit.signatures.base import AbstractSignatureArray
from gambit.util.misc import chunk_slices
from gambit.util.progress import get_progress


def jaccard_generic(set1: Iterable, set2: Iterable) -> float:
	"""Get the Jaccard index of of two arbitrary sets.

	This is primarily used as a slow, pure-Python alternative to :func:`.jaccard_sparse` to be used
	for testing, but can also be used as a generic way to calculate the Jaccard index which works
	with any collection or element type.

	See Also
	--------
	.jaccard_sparse
	.jaccard_bits
	"""
	if not isinstance(set1, Set):
		set1 = set(set1)

	n1 = len(set1)
	n2 = len(set2)
	intersection = len(set1.intersection(set2))
	union = n1 + n2 - intersection

	return 0 if union == 0 else intersection / union


def jaccard_bits(bits1: np.ndarray, bits2: np.ndarray) -> float:
	"""Calculate the Jaccard index between two sets represented as bit arrays ("dense" format for k-mer sets).

	See Also
	--------
	.jaccard_sparse
	"""
	n1 = np.count_nonzero(bits1)
	n2 = np.count_nonzero(bits2)
	intersection = np.count_nonzero(bits1 & bits2)
	union = n1 + n2 - intersection
	return 0 if union == 0 else intersection / union


def jaccard_sparse_array(query: KmerSignature, refs: SignatureArray, out: np.ndarray = None, distance: bool = False) -> np.ndarray:
	"""
	Calculate Jaccard scores between a query k-mer signature and an array of reference signatures in
	``SignatureArray`` format.

	This internally uses Cython code that runs in parallel over all signatures in ``sigarray``.
	Because of Cython limitations ``sigarray.bounds.dtype`` must be ``np.intp``, which is usually
	a 64-bit signed integer. If it is not it will be converted automatically.

	Parameters
	----------
	query : numpy.ndarray
		Query k-mer signature in sparse coordinate format (sorted array of k-mer indices).
	refs : gambit.signatures.SignatureArray
		Array of reference signatures.
	out : Optional[numpy.ndarray]
		Optional pre-allocated array to write results to. Should be the same length as ``sigarray``
		with dtype ``np.float32``.
	distance : bool
		Return Jaccard distances instead of scores.

	Returns
	-------
	numpy.ndarray
		Jaccard score for ``query`` against each element of ``refs``.

	See Also
	--------
	.jaccard_sparse
	.jaccarddist_sparse
	"""
	if out is None:
		out = np.empty(len(refs), SCORE_DTYPE)
	elif out.shape != (len(refs),):
		raise ValueError('Output array length must match signature array.')
	elif out.dtype != SCORE_DTYPE:
		raise ValueError(f'Output array dtype must be {SCORE_DTYPE}, got {out.dtype}')

	values = refs.values
	bounds = refs.bounds.astype(BOUNDS_DTYPE, copy=False)

	_jaccard_sparse_parallel(query, values, bounds, out)

	if distance:
		np.subtract(1, out, out=out)

	return out


def jaccard_sparse_matrix(queries: Sequence[KmerSignature],
                          refs: AbstractSignatureArray,
                          ref_indices: Optional[Sequence[int]] = None,
                          out: Optional[np.ndarray] = None,
                          distance: bool = False,
                          chunksize: Optional[int] = None,
                          progress = None,
                          ) -> np.ndarray:
	"""
	Calculate a Jaccard similarity/distance matrix between an array of query signatures and an
	array of reference signatures.

	The main purpose of this function is to improve querying performance when the reference
	signatures are stored in a file (e.g. using :class:`gambit.signatures.hdf5.HDF5Signatures`)
	by loading them in chunks (via the ``chunksize`` parameter) instead of all in one go.

	Parameters
	----------
	queries
		Query signatures in sparse coordinate format. May be any sequence type, e.g. ``list``.
	refs
		Reference signatures in sparse coordinate format. Must be a type which yields
		``SignatureArray``\\ s when sliced.
	ref_indices
		Optional, indices of ``refs`` to use.
	out
		(Optional) pre-allocated array to write output to.
	distance
		Output Jaccard distances instead of similarities.
	chunksize
		Divide ``refs`` into chunks of this size.
	progress
		Display a progress meter of the number of elements of the output array calculated so far.
		See :func:`gambit.util.progress.get_progress` for a description of allowed values.

	Returns
	-------
	np.ndarray
		Matrix of similarities/distances between query signatures in rows and reference signatures
		in columns.
	"""
	nqueries = len(queries)
	nrefs = len(refs) if ref_indices is None else len(ref_indices)

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
			assert isinstance(ref_chunk, SignatureArray)

			for (i, query) in enumerate(queries):
				jaccard_sparse_array(query, ref_chunk, out=out[i, ref_slice], distance=distance)
				meter.increment(len(ref_chunk))

	return out
