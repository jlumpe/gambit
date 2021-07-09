"""Shared typedefs."""

cimport numpy as np


# Type for similarity scores
ctypedef np.float32_t SCORE_T

# Type for bounds on c_jaccard_coords_col
ctypedef np.intp_t BOUNDS_T

# Fused type for storing k-mer coordinates/indices
ctypedef fused COORDS_T:
	np.int16_t
	np.uint16_t
	np.int32_t
	np.uint32_t
	np.int64_t
	np.uint64_t

# Copy of COORDS_T, used when two arguments have types in this set but may be different than each other.
ctypedef fused COORDS_T_2:
	np.int16_t
	np.uint16_t
	np.int32_t
	np.uint32_t
	np.int64_t
	np.uint64_t
