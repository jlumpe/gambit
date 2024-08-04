"""Cython module for working with DNA sequences and k-mers.

Note: each of the 4 Python functions here have a C counterpart that does the actual work. The Python
version is just a wrapper that does any needed conversion, allocates buffers, and raises exceptions
if needed. The separation currently isn't necessary as the C functions aren't used anywhere else
outside the wrappers, but they may be in the future. Handling exceptions in the Python wrappers only
allows the C functions to be declared with nogil.
"""


def kmer_to_index(const CHAR[:] kmer):
	"""kmer_to_index(kmer: bytes) -> int

	Convert k-mer byte string to its integer index.
	"""
	cdef:
		uint64_t idx
		bint exc = False

	if kmer.shape[0] > 32:
		raise ValueError('k must be <= 32')

	idx = c_kmer_to_index(kmer, &exc)

	if exc:
		raise ValueError('Invalid character in k-mer')

	return idx


cdef uint64_t c_kmer_to_index(const CHAR[:] kmer, bint *exc) nogil:
	cdef:
		uint64_t idx = 0
		int i, k = kmer.shape[0]
		CHAR nuc

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
			exc[0] = True
			return 0

	return idx


def kmer_to_index_rc(const CHAR[:] kmer):
	"""kmer_to_index_rc(kmer: bytes) -> int

	Get the integer index of the reverse complement of a k-mer byte string.
	"""
	cdef:
		uint64_t idx
		bint exc = False

	if kmer.shape[0] > 32:
		raise ValueError('k must be <= 32')

	idx = c_kmer_to_index_rc(kmer, &exc)

	if exc:
		raise ValueError('Invalid character in k-mer')

	return idx


cdef uint64_t c_kmer_to_index_rc(const CHAR[:] kmer, bint *exc) nogil:
	cdef:
		uint64_t idx = 0
		int i, k = kmer.shape[0]
		CHAR nuc

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
			exc[0] = True
			return 0

	return idx


def index_to_kmer(index, int k):
	"""index_to_kmer(index: int, kmer: int) -> bytes

	Convert k-mer index to sequence.
	"""
	buf = bytearray(k)
	c_index_to_kmer(index, buf)
	return bytes(buf)


cdef void c_index_to_kmer(uint64_t index, CHAR[:] out) nogil:
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
