"""Cython module for working with DNA sequences and k-mers."""

cimport numpy as np

from .types cimport COORDS_T

ctypedef unsigned char CHAR


cdef np.uint64_t c_kmer_to_index(const CHAR[:]) nogil except? 0
cdef np.uint64_t c_kmer_to_index_rc(const CHAR[:]) nogil except? 0
cdef void c_index_to_kmer(COORDS_T, int, char*) nogil
cdef inline char nuc_complement(char) nogil
cdef void c_revcomp(const char*, int, char*) nogil
