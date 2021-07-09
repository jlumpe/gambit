"""Cython functions for calculating k-mer similarity metrics."""

from .types cimport SCORE_T, BOUNDS_T, COORDS_T, COORDS_T_2


cdef SCORE_T c_jaccard_sparse(COORDS_T[:] coords1, COORDS_T_2[:] coords2) nogil
