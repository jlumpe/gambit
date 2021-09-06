"""Core functions for searching for and working with k-mers.

Note that all code in this module operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

from typing import Sequence, Optional, Union, NewType, Dict, Any

import numpy as np
from attr import attrs, attrib

from gambit._cython.kmers import kmer_to_index, index_to_kmer, reverse_complement
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


def coords_dtype(k : int) -> np.dtype:
	"""Get the smallest unsigned integer dtype that can store k-mer indices for the given ``k``.

	Parameters
	----------
	k : int

	Returns
	-------
	numpy.dtype
	"""
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
	coords_dtype
		Smallest unsigned integer dtype that can store k-mer indices.
	"""
	k: int = attrib()
	prefix: bytes = attrib(
		converter=lambda v: v.upper().encode('ascii') if isinstance(v, str) else v,
	)
	prefix_str: str
	prefix_len: int
	total_len: int
	idx_len: int
	coords_dtype: np.dtype

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
		object.__setattr__(self, 'idx_len', 4 ** self.k)
		object.__setattr__(self, 'coords_dtype', coords_dtype(self.k))

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


def find_kmers(kspec: KmerSpec,
               seq: Union[bytes, str],
               *,
               sparse: bool = True,
               dense_out: Optional[np.ndarray] = None,
               ) -> KmerSignature:
	"""Find all k-mers in a DNA sequence.

	Searches sequence both backwards and forwards (reverse complement). The sequence may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kspec : .KmerSpec
		K-mer spec to use for search.
	seq
		Sequence to search within as ``bytes`` or ``str``. If ``str`` will be encoded as ASCII.
		Lower-case characters are OK and will be matched as upper-case.
	dense_out : numpy.ndarray
		Pre-allocated numpy array to write dense output to. Should be of length ``kspec.idx_len``.
		Note that this is still used as working space even if ``sparse=True``. Should be zeroed
		prior to use (although if not the result will effectively be the bitwise AND between its
		previous value and k-mers found in ``data``.
	sparse : bool
		If True return k-mers in sparse coordinate format rather than dense (bit vector) format.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	gambit.io.seq.find_kmers_parse
	"""
	if dense_out is None:
		dense_out = np.zeros(kspec.idx_len, dtype=bool)

	# Convert sequence to bytes
	if not isinstance(seq, bytes):
		if not isinstance(seq, str):
			seq = str(seq)

		seq = seq.encode('ascii')

	# Convert to upper-case only if needed
	nucs_lower = set(NUCLEOTIDES.lower())
	for char in seq:
		if char in nucs_lower:
			seq = seq.upper()
			break

	_find_kmers(kspec, seq, dense_out)

	if sparse:
		return dense_to_sparse(dense_out)
	else:
		return dense_out


def _find_kmers(kspec, seq, out):
	"""Actual implementation of find_kmers.

	Parameters
	----------
	kspec : KmerSpec
	seq : bytes
		Upper-case ASCII nucleotide codes.
	out : np.ndarray
		Write dense output to this array.
	"""

	# Reverse complement of prefix
	rcprefix = reverse_complement(kspec.prefix)

	# Search forward
	start = 0
	while True:
		loc = seq.find(kspec.prefix, start, -kspec.k)
		if loc < 0:
			break

		kmer = seq[loc + kspec.prefix_len:loc + kspec.total_len]
		if not isinstance(kmer, bytes):
			kmer = str(kmer).encode('ascii')

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1

	# Search backward
	start = kspec.k
	while True:
		loc = seq.find(rcprefix, start)
		if loc < 0:
			break

		rckmer = seq[loc - kspec.k:loc]
		if not isinstance(rckmer, bytes):
			rckmer = str(rckmer).encode('ascii')
		kmer = reverse_complement(rckmer)

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1


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
	idx_len = k_or_kspec.idx_len if isinstance(k_or_kspec, KmerSpec) else 4 ** k_or_kspec
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

	range_ = 4 ** (from_kspec.k - extra_len)
	start = extra_ind * range_
	stop = (extra_ind + 1) * range_
	reduce = from_kspec.k - to_kspec.k - extra_len

	return start, stop, reduce


def convert_dense(from_kspec: KmerSpec, to_kspec: KmerSpec, vec: np.ndarray) -> np.ndarray:
	"""Convert a k-mer signature in dense format from one ``KmerSpec`` to another.

	In the ideal case, if ``vec`` is the result of ``find_kmers(from_kspec, seq, sparse=False)``
	the output of this function should be identical to ``find_kmers(to_kspec, seq, sparse=False)``\\ .
	In reality this may not hold if any potential matches of ``from_kspec`` in ``seq`` are discarded
	due to an invalid nucleotide which is not included in the corresponding ``to_kspec`` match.
	"""
	check_can_convert(from_kspec, to_kspec)
	start, stop, reduce = _convert_params(from_kspec, to_kspec)
	block_size = 4 ** reduce

	out = np.zeros(to_kspec.idx_len, dtype=bool)

	for i in range(block_size):
		out |= vec[start+i:stop:block_size]

	return out


def convert_sparse(from_kspec: KmerSpec, to_kspec: KmerSpec, sig: KmerSignature) -> KmerSignature:
	"""Convert a k-mer signature in sparse format from one ``KmerSpec`` to another.

	In the ideal case, if ``sig`` is the result of ``find_kmers(from_kspec, seq)``
	the output of this function should be identical to ``find_kmers(to_kspec, seq)``\\ .
	In reality this may not hold if any potential matches of ``from_kspec`` in ``seq`` are discarded
	due to an invalid nucleotide which is not included in the corresponding ``to_kspec`` match.
	"""
	assert can_convert(from_kspec, to_kspec)
	start, stop, reduce = _convert_params(from_kspec, to_kspec)
	reduce_bits = 2 * reduce

	out = np.empty(len(sig), dtype=to_kspec.coords_dtype)
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
