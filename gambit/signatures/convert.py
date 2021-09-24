"""Convert signatures between representations or from one KmerSpec to another."""

from typing import Sequence, Union

import numpy as np

from .base import KmerSignature
from gambit.kmers import KmerSpec, nkmers, kmer_to_index


def dense_to_sparse(vec: Sequence[bool]) -> KmerSignature:
	"""Convert k-mer set from dense bit vector to sparse coordinate representation.

	Parameters
	----------
	vec : numpy.ndarray
		Boolean vector indicating which k-mers are present.

	Returns
	-------
	numpy.ndarray
		Sorted array  of coordinates of k-mers present in vector. Data type will be ``numpy.intp``.

	See Also
	--------
	.sparse_to_dense
	"""
	return np.flatnonzero(vec)


def sparse_to_dense(k_or_kspec: Union[int, KmerSpec],  coords: KmerSignature) -> np.ndarray:
	"""Convert k-mer set from sparse coordinate representation back to dense bit vector.

	Parameters
	----------
	k_or_kspec : Union[int, KmerSpec]
		Value of k or a :class:`.KmerSpec` instance.
	coords : numpy.ndarray
		Sparse coordinate array.

	Returns
	-------
	numpy.ndarray
		Dense k-mer bit vector.

	See Also
	--------
	.dense_to_sparse
	"""
	idx_len = k_or_kspec.nkmers if isinstance(k_or_kspec, KmerSpec) else nkmers(k_or_kspec)
	vec = np.zeros(idx_len, dtype=np.bool_)
	vec[coords] = 1
	return vec


def can_convert(from_kspec: KmerSpec, to_kspec: KmerSpec) -> bool:
	"""Check if signatures from one KmerSpec can be converted to another.

	Conversion is possible if ``to_kspec.prefix`` is equal to or starts with ``from_kspec.prefix``
	and ``to_kspec.total_len <= from_kspec.total_len``.
	"""
	return to_kspec.prefix.startswith(from_kspec.prefix) and to_kspec.total_len <= from_kspec.total_len


def check_can_convert(from_kspec: KmerSpec, to_kspec: KmerSpec):
	"""
	Check that signatures can be converted from one KmerSpec to another or raise an error with an
	informative message.

	Raises
	------
	ValueError
		If conversion is not possible.
	"""
	if not to_kspec.prefix.startswith(from_kspec.prefix):
		raise ValueError('Destination prefix must start with source prefix.')
	if to_kspec.total_len > from_kspec.total_len:
		raise ValueError('Cannot convert to KmerSpec with longer total length.')


def _convert_params(from_kspec: KmerSpec, to_kspec: KmerSpec):
	extra_prefix = to_kspec.prefix[from_kspec.prefix_len:]
	extra_ind = kmer_to_index(extra_prefix)
	extra_len = len(extra_prefix)

	range_ = nkmers(from_kspec.k - extra_len)
	start = extra_ind * range_
	stop = (extra_ind + 1) * range_
	reduce = from_kspec.k - to_kspec.k - extra_len

	return start, stop, reduce


def convert_dense(from_kspec: KmerSpec, to_kspec: KmerSpec, vec: np.ndarray) -> np.ndarray:
	"""Convert a k-mer signature in dense format from one ``KmerSpec`` to another.

	In the ideal case, if ``vec`` is the result of ``calc_signature(from_kspec, seq, sparse=False)``
	the output of this function should be identical to ``calc_signature(to_kspec, seq, sparse=False)``.
	In reality this may not hold if any potential matches of ``from_kspec`` in ``seq`` are discarded
	due to an invalid nucleotide which is not included in the corresponding ``to_kspec`` match.
	"""
	check_can_convert(from_kspec, to_kspec)
	start, stop, reduce = _convert_params(from_kspec, to_kspec)
	block_size = nkmers(reduce)

	out = np.zeros(to_kspec.nkmers, dtype=bool)

	for i in range(block_size):
		out |= vec[start+i:stop:block_size]

	return out


def convert_sparse(from_kspec: KmerSpec, to_kspec: KmerSpec, sig: KmerSignature) -> KmerSignature:
	"""Convert a k-mer signature in sparse format from one ``KmerSpec`` to another.

	In the ideal case, if ``sig`` is the result of ``calc_signature(from_kspec, seq)``
	the output of this function should be identical to ``calc_signature(to_kspec, seq)``.
	In reality this may not hold if any potential matches of ``from_kspec`` in ``seq`` are discarded
	due to an invalid nucleotide which is not included in the corresponding ``to_kspec`` match.
	"""
	assert can_convert(from_kspec, to_kspec)
	start, stop, reduce = _convert_params(from_kspec, to_kspec)
	reduce_bits = 2 * reduce

	out = np.empty(len(sig), dtype=to_kspec.index_dtype)
	i = 0
	next_ = start

	for from_idx in sig:
		if from_idx < next_:
			continue
		if from_idx >= stop:
			break

		to_idx = (from_idx - start) >> reduce_bits
		out[i] = to_idx
		i += 1

		# Next possible input index that won't reduce to the same output
		next_ = ((to_idx + 1) << reduce_bits) + start

	out.resize(i)
	return out