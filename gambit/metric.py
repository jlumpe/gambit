"""Calculate the Jaccard index/distance between sets."""

from collections.abc import Set
from typing import Iterable

import numpy as np

from gambit._cython.metric import BOUNDS_DTYPE, SCORE_DTYPE, jaccard_sparse, jaccarddist_sparse, \
	_jaccard_sparse_parallel
from gambit.kmers import KmerSignature
from gambit.signatures import SignatureArray


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
