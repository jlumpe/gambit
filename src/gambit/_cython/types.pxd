"""Shared typedefs."""

from libc.stdint cimport uint16_t, uint32_t, uint64_t, intptr_t


# Type for similarity scores
ctypedef float SCORE_T

# Type for bounds on c_jaccard_coords_col
# This should be equal to Numpy's intp dtype
ctypedef intptr_t BOUNDS_T

# Fused type for storing k-mer coordinates/indices
ctypedef fused COORDS_T:
	uint16_t
	uint32_t
	uint64_t

# Copy of COORDS_T, used when two arguments have types in this set but may be different than each other.
ctypedef fused COORDS_T_2:
	uint16_t
	uint32_t
	uint64_t
