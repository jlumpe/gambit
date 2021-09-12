"""Core functions for searching for and working with k-mers.

Note that all code in this module operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

from typing import Sequence, Union, NewType, Dict, Any

import numpy as np
from attr import attrs, attrib

from gambit._cython.kmers import kmer_to_index, kmer_to_index_rc, index_to_kmer, revcomp
from gambit.io.json import Jsonable


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'


#: Type for k-mer signatures (k-mer sets in sparse coordinate format)
KmerSignature = NewType('KmerSignature', np.ndarray)
# TODO - use nptyping package to specify dimensions and data type?


def validate_dna_seq_bytes(seq : bytes):
	"""Check that a sequence contains only valid nucleotide codes.

	Parameters
	----------
	seq : bytes
		ASCII-encoded nucleotide sequence.

	Raises
	------
	ValueError
		If the sequence contains an invalid nucleotide.
	"""
	for i, nuc in enumerate(seq):
		if nuc not in NUCLEOTIDES:
			raise ValueError(f'Invalid byte at position {i}: {nuc}')


def nkmers(k: int) -> int:
	"""Get the number of possible distinct k-mers for a given value of ``k``."""
	return 4 ** k


def index_dtype(k: int) -> np.dtype:
	"""Get the smallest unsigned integer dtype that can store k-mer indices for the given ``k``."""
	if k <= 4:
		return np.dtype('u1')
	elif k <= 8:
		return np.dtype('u2')
	elif k <= 16:
		return np.dtype('u4')
	elif k <= 32:
		return np.dtype('u8')
	else:
		return None


@attrs(frozen=True, repr=False, cmp=False)
class KmerSpec(Jsonable):
	"""Specifications for a k-mer search operation.

	Parameters
	----------
	k : int
		Value of :attr:`k` attribute.
	prefix : str or bytes
		Value of :attr:`prefix` attribute. If ``str`` and not ``bytes`` will be encoded as ascii.

	Attributes
	----------
	k
		Number of nucleotides in k-mer *after* prefix.
	prefix
		Constant prefix of k-mers to search for, upper-case nucleotide codes
		as ascii-encoded ``bytes``.
	prefix_str
		Prefix as string.
	prefix_len
		Number of nucleotides in prefix.
	total_len
		Sum of ``prefix_len`` and ``k``.
	idx_len
		Maximum value (plus one) of integer needed to index one of the
		found k-mers. Also the number of possible k-mers fitting the spec.
		Equal to ``4 ** k``.
	index_dtype
		Smallest unsigned integer dtype that can store k-mer indices.
	"""
	k: int = attrib()
	prefix: bytes = attrib(
		converter=lambda v: v.upper().encode('ascii') if isinstance(v, str) else v,
	)
	prefix_str: str
	prefix_len: int
	total_len: int
	nkmers: int
	index_dtype: np.dtype

	@k.validator
	def _validate_k(self, attribute, value):
		if value < 1:
			raise ValueError('k must be positive')

	@prefix.validator
	def _validate_prefix(self, attribute, value):
		validate_dna_seq_bytes(value)

	def __attrs_post_init__(self):
		object.__setattr__(self, 'prefix_str', self.prefix.decode('ascii'))
		object.__setattr__(self, 'prefix_len', len(self.prefix))
		object.__setattr__(self, 'total_len', self.k + self.prefix_len)
		object.__setattr__(self, 'nkmers', nkmers(self.k))
		object.__setattr__(self, 'index_dtype', index_dtype(self.k))

	def __get_newargs__(self):
		return self.k, self.prefix

	def __eq__(self, other):
		return isinstance(other, KmerSpec) and\
			self.k == other.k and\
			self.prefix == other.prefix

	def __repr__(self):
		return f'{type(self).__name__}({self.k}, {self.prefix_str!r})'

	def __to_json__(self):
		return dict(k=int(self.k), prefix=self.prefix_str)

	@classmethod
	def __from_json__(cls, data: Dict[str, Any]) -> 'KmerSpec':
		return cls(data['k'], data['prefix'])


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
	and ``to_kspec.total_len <= from_kspec.total_len``\\ .
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
	the output of this function should be identical to ``calc_signature(to_kspec, seq, sparse=False)``\\ .
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
	the output of this function should be identical to ``calc_signature(to_kspec, seq)``\\ .
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
