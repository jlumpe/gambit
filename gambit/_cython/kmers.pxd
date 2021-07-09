"""Cython module for working with DNA sequences and k-mers."""

cimport numpy as np

from .types cimport COORDS_T


cdef np.uint32_t c_kmer_to_index32(const char*, int) except? 0
cdef np.uint64_t c_kmer_to_index64(const char*, int) except? 0
cdef void c_index_to_kmer(COORDS_T, int, char*) nogil
cdef inline char nuc_complement(char) nogil
cdef void c_reverse_complement(const char*, int, char*) nogil
