from .types cimport SCORE_T, BOUNDS_T, COORDS_T, COORDS_T_2

cdef SCORE_T c_jaccarddist(COORDS_T[:] coords1, COORDS_T_2[:] coords2) nogil
