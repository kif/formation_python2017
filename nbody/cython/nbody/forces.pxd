cimport numpy as cnp

cdef cnp.float64_t[::1] force(cnp.float64_t[::1] p1,
                              cnp.float64_t[::1] p2,
                              cnp.float64_t m2)
