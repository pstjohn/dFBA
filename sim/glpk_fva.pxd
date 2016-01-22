import numpy as np
cimport numpy as np

from dFBA.sim.glpk cimport GLPKfba

cdef class GLPKfva(GLPKfba):

    cdef: 
        np.float_t[:] c_inner
        np.float_t[:] c_outer
        np.float_t[:] lb_orig
        np.float_t[:] ub_orig

        double fva_tol

    cpdef int solve(self)


