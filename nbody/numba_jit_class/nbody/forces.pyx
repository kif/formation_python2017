from libc.math cimport sqrt
from .physics import gamma_si as GAMMA_SI, eps as EPS
cimport numpy as cnp
import numpy as np

cdef: #just a couple of constants
    cnp.float64_t gamma_si=GAMMA_SI
    cnp.float64_t eps = EPS

cdef cnp.float64_t[::1] force(cnp.float64_t[::1] p1,
                              cnp.float64_t[::1] p2,
                              cnp.float64_t m2):
    cdef:
        cnp.float64_t dx, dy, dist, F
        cnp.float64_t[::1] res = np.zeros(2, dtype=np.float64)
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dist = sqrt(dx**2 + dy**2 + eps)
    if dist > 0:
        F = (gamma_si * m2) / (dist*dist*dist)
        res[0] = F * dx
        res[1] = F * dy
    return res
