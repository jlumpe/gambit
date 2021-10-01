"""Core functions for searching for and working with k-mers."""

from typing import Dict, Any, Iterator

import numpy as np
from attr import attrs, attrib

import gambit._cython.kmers as ckmers
from gambit._cython.kmers import index_to_kmer
from gambit.seq import NUCLEOTIDES, DNASeq, seq_to_bytes, validate_dna_seq_bytes, revcomp
from gambit.io.json import Jsonable


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


def kmer_to_index(kmer: DNASeq) -> int:
	"""Convert a k-mer to its integer index.

	Raises
	------
	ValueError
		If an invalid nucleotide code is encountered.
	"""
	return ckmers.kmer_to_index(seq_to_bytes(kmer))


def kmer_to_index_rc(kmer: DNASeq) -> int:
	"""Get the integer index of a k-mer's reverse complement.

	Raises
	------
	ValueError
		If an invalid nucleotide code is encountered.
	"""
	return ckmers.kmer_to_index_rc(seq_to_bytes(kmer))


@attrs(frozen=True, repr=False, init=False)
class KmerSpec(Jsonable):
	"""Specifications for a k-mer search operation.

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
	prefix: bytes = attrib()
	prefix_str: str = attrib(eq=False)
	prefix_len: int = attrib(eq=False)
	total_len: int = attrib(eq=False)
	nkmers: int = attrib(eq=False)
	index_dtype: np.dtype = attrib(eq=False)

	def __init__(self, k: int, prefix: DNASeq):
		"""
		Parameters
		----------
		k
			Value of :attr:`k` attribute.
		prefix
			Value of :attr:`prefix` attribute. Will be converted to ``bytes``.
		"""
		if k < 1:
			raise ValueError('k must be positive')

		prefix = seq_to_bytes(prefix).upper()
		validate_dna_seq_bytes(prefix)

		self.__attrs_init__(
			k=k,
			prefix=prefix,
			prefix_str=prefix.decode('ascii'),
			prefix_len=len(prefix),
			total_len=k + len(prefix),
			nkmers=nkmers(k),
			index_dtype=index_dtype(k),
		)

	def __get_newargs__(self):
		return self.k, self.prefix

	def __repr__(self):
		return f'{type(self).__name__}({self.k}, {self.prefix_str!r})'

	def __to_json__(self):
		return dict(k=int(self.k), prefix=self.prefix_str)

	@classmethod
	def __from_json__(cls, data: Dict[str, Any]) -> 'KmerSpec':
		return cls(data['k'], data['prefix'])


@attrs(slots=True)
class KmerMatch:
	"""Represents a

	Attributes
	----------
	kmerspec
		K-mer spec used for search.
	seq
		The sequence searched within.
	pos
		Index of first nucleotide of prefix in ``seq``.
	reverse
		If the match is on the reverse strand.
	"""
	kmerspec: KmerSpec = attrib()
	seq: DNASeq = attrib()
	pos: int = attrib()
	reverse: bool = attrib()

	def kmer_indices(self) -> slice:
		"""Index range for k-mer in sequence (without prefix)."""
		if self.reverse:
			return slice(self.pos - self.kmerspec.total_len + 1, self.pos - self.kmerspec.prefix_len + 1)
		else:
			return slice(self.pos + self.kmerspec.prefix_len, self.pos + self.kmerspec.total_len)

	def full_indices(self) -> slice:
		"""Index range for prefix plus k-mer in sequence."""
		if self.reverse:
			return slice(self.pos - self.kmerspec.total_len + 1, self.pos + 1)
		else:
			return slice(self.pos, self.pos + self.kmerspec.total_len)

	def kmer(self) -> bytes:
		"""Get matched k-mer sequence."""
		kmer = seq_to_bytes(self.seq[self.kmer_indices()])
		return revcomp(kmer) if self.reverse else kmer

	def kmer_index(self) -> int:
		"""Get index of matched k-mer.

		Raises
		------
		ValueError
			If the k-mer contains invalid nucleotides.
		"""
		kmer = self.seq[self.kmer_indices()]
		return kmer_to_index_rc(kmer) if self.reverse else kmer_to_index(kmer)


def find_kmers(kmerspec: KmerSpec, seq: DNASeq) -> Iterator[KmerMatch]:
	"""Locate k-mers with the given prefix in a DNA sequence.

	Searches sequence both backwards and forwards (reverse complement). The sequence may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kmerspec
		K-mer spec to use for search.
	seq
		Sequence to search within. Lowercase characters are OK and will be matched as uppercase.

	Returns
	-------
	Iterator[KmerMatch]
		Iterator of :class:`.KmerMatch` objects.
	"""

	haystack = seq_to_bytes(seq)

	# Convert to uppercase only if needed
	nucs_lower = NUCLEOTIDES.lower()
	for char in haystack:
		if char in nucs_lower:
			haystack = haystack.upper()
			break

	# Find forward
	start = 0

	while True:
		loc = haystack.find(kmerspec.prefix, start, -kmerspec.k)
		if loc < 0:
			break

		yield KmerMatch(kmerspec, seq, loc, False)

		start = loc + 1

	# Find reverse
	prefix_rc = revcomp(kmerspec.prefix)
	start = kmerspec.k

	while True:
		loc = haystack.find(prefix_rc, start)
		if loc < 0:
			break

		yield KmerMatch(kmerspec, seq, loc + kmerspec.prefix_len - 1, True)

		start = loc + 1
