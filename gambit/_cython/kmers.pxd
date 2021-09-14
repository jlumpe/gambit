# cython: language_level = 3str

cimport numpy as np

ctypedef unsigned char CHAR


cpdef np.uint64_t kmer_to_index(const CHAR[:]) nogil except? 0
cpdef np.uint64_t kmer_to_index_rc(const CHAR[:]) nogil except? 0
cdef void c_index_to_kmer(np.uint64_t, CHAR[:]) nogil
cdef void c_revcomp(const CHAR[:], CHAR[:]) nogil
