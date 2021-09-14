# cython: language_level = 3str, wraparound = False, boundscheck = False

"""Cython module for working with DNA sequences and k-mers."""


cpdef np.uint64_t kmer_to_index(const CHAR[:] kmer) nogil except? 0:
	"""kmer_to_index(kmer)

	Convert k-mer byte string to its integer index.
	"""
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


cpdef np.uint64_t kmer_to_index_rc(const CHAR[:] kmer) nogil except? 0:
	"""kmer_to_index_rc(kmer)

	Get the integer index of the reverse complement of a k-mer byte string.
	"""
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


def index_to_kmer(index, int k):
	"""index_to_kmer(index: int, kmer: int) -> bytes

	Convert k-mer index to sequence.
	"""
	buf = bytearray(k)
	c_index_to_kmer(index, buf)
	return bytes(buf)


cdef void c_index_to_kmer(np.uint64_t index, CHAR[:] out) nogil:
	"""Convert k-mer index to sequence."""
	cdef:
		int k = out.shape[0]
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


def revcomp(const CHAR[:] seq):
	"""revcomp(seq: bytes) -> bytes

	Get the reverse complement of a nucleotide sequence.

	Parameters
	----------
	seq : bytes
		ASCII-encoded nucleotide sequence. Case does not matter.

	Returns
	-------
	bytes
		Reverse complement sequence. All characters in the input which are not valid nucleotide
		codes will appear unchanged in the corresponding reverse position.
	"""
	buf = bytearray(len(seq))
	c_revcomp(seq, buf)
	return bytes(buf)


cdef void c_revcomp(const CHAR[:] seq, CHAR[:] out) nogil:
	"""Get the reverse complement of a nucleotide sequence.

	Parameters
	----------
	seq
		Sequence as byte array.
	out
		Byte buffer of same length as seq to write output to.
	"""

	cdef:
		int i, n = seq.shape[0]
		CHAR nuc, nuc2

	for i in range(n):
		nuc = seq[i]

		if nuc == 'A':
			nuc2 =  'T'
		elif nuc == 'a':
			nuc2 =  't'
		elif nuc == 'T':
			nuc2 =  'A'
		elif nuc == 't':
			nuc2 =  'a'
		elif nuc == 'G':
			nuc2 =  'C'
		elif nuc == 'g':
			nuc2 =  'c'
		elif nuc == 'C':
			nuc2 =  'G'
		elif nuc == 'c':
			nuc2 =  'g'
		else:
			nuc2 =  nuc

		out[n - i - 1] = nuc2
