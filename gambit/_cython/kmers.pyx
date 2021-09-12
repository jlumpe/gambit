# cython: language_level = 3str, wraparound = False

"""Cython module for working with DNA sequences and k-mers."""

import numpy as np

from libc.stdlib cimport malloc, free


# Allowable Numpy dtypes for coordinate arrays
COORDS_DTYPES = frozenset(map(np.dtype, [
	'i2',
	'u2',
	'i4',
	'u4',
	'i8',
	'u8',
]))


def kmer_to_index(kmer):
	"""kmer_to_index(kmer)

	Convert k-mer byte string to its integer index.

	Parameters
	----------
	kmer : Union[str, bytes, bytearray]
		K-mer as string or bytes.

	Returns
	-------
	int
		Index of k-mer.

	Raises
	------
	ValueError
		If an invalid nucleotide code is encountered.
	"""
	if isinstance(kmer, str):
		kmer = kmer.encode('ascii')
	if isinstance(kmer, (bytes, bytearray)):
		return c_kmer_to_index(kmer)
	raise TypeError(f'Expected str, bytes, or bytearray, got {type(kmer)}')

cdef np.uint64_t c_kmer_to_index(const CHAR[:] kmer) nogil except? 0:
	cdef:
		np.uint64_t idx = 0
		int i, k = kmer.shape[0]
		CHAR nuc

	if k > 32:
		raise ValueError('k must be <= 32')

	for i in range(k):
		nuc = kmer[i]

		idx <<= 2

		nuc &= 0b11011111  # To upper case
		if nuc == 'A':
			idx += 0
		elif nuc == 'C':
			idx += 1
		elif nuc == 'G':
			idx += 2
		elif nuc == 'T':
			idx += 3
		else:
			raise ValueError(nuc)

	return idx


def kmer_to_index_rc(kmer):
	"""kmer_to_index_rc(kmer)

	Get the integer index of the reverse complement of a k-mer.

	Parameters
	----------
	kmer : Union[str, bytes, bytearray]
		K-mer as string or bytes.

	Returns
	-------
	int
		Index of k-mer's reverse complement.

	Raises
	------
	ValueError
		If an invalid nucleotide code is encountered.
	"""
	if isinstance(kmer, str):
		kmer = kmer.encode('ascii')
	if isinstance(kmer, (bytes, bytearray)):
		return c_kmer_to_index_rc(kmer)
	raise TypeError(f'Expected str, bytes, or bytearray, got {type(kmer)}')

cdef np.uint64_t c_kmer_to_index_rc(const CHAR[:] kmer) nogil except? 0:
	cdef:
		np.uint64_t idx = 0
		int i, k = kmer.shape[0]
		CHAR nuc

	if k > 32:
		raise ValueError('k must be <= 32')

	for i in range(k):
		nuc = kmer[k - i - 1]

		idx <<= 2

		nuc &= 0b11011111  # To upper case
		if nuc == 'A':
			idx += 3
		elif nuc == 'C':
			idx += 2
		elif nuc == 'G':
			idx += 1
		elif nuc == 'T':
			idx += 0
		else:
			raise ValueError(nuc)

	return idx


def index_to_kmer(np.uint64_t index, int k):
	"""index_to_kmer(index, k)

	Convert k-mer index to sequence.

	Parameters
	----------
	index : int
		K-mer index.
	k : int
		Length of k-mer.

	Returns
	-------
	bytes
		K-mer as byte string.
	"""

	cdef char* buf = <char*>malloc((k + 1) * sizeof(char))

	try:

		c_index_to_kmer(index, k, buf)
		buf[k] = 0  # Null-terminate it
		return <bytes>buf

	finally:
		free(buf)


def reverse_complement(bytes seq):
	"""reverse_complement(seq)

	Get the reverse complement of a nucleotide sequence.

	Parameters
	----------
	seq : bytes
		ASCII-encoded nucleotide sequence. Case does not matter.

	Returns
	-------
	bytes
		Reverse complement sequence. All characters in the input which are not valid nucleotide codes will appear
		unchanged in the corresponding reverse position.
	"""

	cdef:
		int l = len(seq)
		char* buf = <char*>malloc((l + 1) * sizeof(char))

	try:
		c_reverse_complement(seq, l, buf)
		buf[l] = 0  # Null-terminate it
		return <bytes>buf

	finally:
		free(buf)


cdef void c_index_to_kmer(COORDS_T index, int k, char* out) nogil:
	"""Convert k-mer index to sequence.

	Output will be in upper-case characters.

	Parameters
	----------
	index
		Index of k-mer
	k
		Length of k-mer
	out
		Destination buffer of length k sequence will be written to.
	"""
	cdef:
		int i, nuc_index
		char nuc

	for i in range(k):

		nuc_index = index % 4

		if nuc_index == 0:
			nuc = 'A'
		elif nuc_index == 1:
			nuc = 'C'
		elif nuc_index == 2:
			nuc = 'G'
		else:
			nuc = 'T'

		out[k - i - 1] = nuc

		index >>= 2


cdef inline char nuc_complement(char nuc) nogil:
	"""Get the complement of a nucleotide.

	If the input is a valid nucleotide code the output will be the code of the
	complement nucleotide in the same case. If the input is not a valid
	nucleotide code it will be returned unchanged.
	"""
	if nuc == 'A':
		return 'T'
	elif nuc == 'a':
		return 't'
	elif nuc == 'T':
		return 'A'
	elif nuc == 't':
		return 'a'
	elif nuc == 'G':
		return 'C'
	elif nuc == 'g':
		return 'c'
	elif nuc == 'C':
		return 'G'
	elif nuc == 'c':
		return 'g'
	else:
		return nuc


cdef void c_reverse_complement(const char* seq, int l, char* out) nogil:
	"""Get the reverse complement of a nucleotide sequence.

	Parameters
	----------
	seq
		Pointer to start of sequence
	l
		Length of sequence
	out
		Pointer to buffer of same length as seq to store the output.
	"""

	cdef int i

	for i in range(l):
		out[l - i - 1] = nuc_complement(seq[i])
