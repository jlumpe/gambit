cimport numpy as np

ctypedef unsigned char CHAR


cdef np.uint64_t c_kmer_to_index(const CHAR[:], bint*) nogil
cdef np.uint64_t c_kmer_to_index_rc(const CHAR[:], bint*) nogil
cdef void c_index_to_kmer(np.uint64_t, CHAR[:]) nogil
cdef void c_revcomp(const CHAR[:], CHAR[:]) nogil
