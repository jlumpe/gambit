"""Core functions for searching for and working with k-mers.

Note that all code in this module operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

from typing import Sequence, Optional, Union, NewType

import numpy as np
from attr import attrs, attrib

from gambit._cython.kmers import kmer_to_index, index_to_kmer, reverse_complement
from gambit.io.json import Jsonable


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'


#: Type for k-mer signatures (k-mer sets in sparse coordinate format)
# TODO - use nptyping package to specify dimensions and data type?
KmerSignature = NewType('KmerSignature', np.ndarray)


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
	prefix : bytes
		Constant prefix of k-mers to search for, upper-case nucleotide codes
		as ascii-encoded ``bytes``.
	k : int
		Number of nucleotides in k-mer *after* prefix.
	prefix_len : int
		Number of nucleotides in prefix.
	total_len : int
		Sum of ``prefix_len`` and ``k``.
	idx_len : int
		Maximum value (plus one) of integer needed to index one of the
		found k-mers. Also the number of possible k-mers fitting the spec.
		Equal to ``4 ** k``.
	coords_dtype : numpy.dtype
		Smallest unsigned integer dtype that can store k-mer indices.
	"""
	k: int = attrib()
	prefix: bytes = attrib(
		converter=lambda v: v.upper().encode('ascii') if isinstance(v, str) else v,
	)
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
		return f'{type(self).__name__}({self.k}, {self.prefix.decode("ascii")!r})'

	def __to_json__(self):
		return dict(k=self.k, prefix=self.prefix.decode('ascii'))

	@classmethod
	def __from_json__(cls, data):
		return cls(data['k'], data['prefix'])


def find_kmers(
               kspec: KmerSpec,
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
		:func:`gambit.kmers.vec_to_coords`).

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
