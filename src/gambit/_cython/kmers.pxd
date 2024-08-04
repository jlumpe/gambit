from libc.stdint cimport uint64_t, intptr_t

ctypedef unsigned char CHAR


cdef uint64_t c_kmer_to_index(const CHAR[:], bint*) nogil
cdef uint64_t c_kmer_to_index_rc(const CHAR[:], bint*) nogil
cdef void c_index_to_kmer(uint64_t, CHAR[:]) nogil
cdef void c_revcomp(const CHAR[:], CHAR[:]) nogil
