import numpy as np

cdef class GLPKfva(GLPKfba):

    def __init__(self, cobra_model, 
                 np.ndarray[np.float_t, ndim=1] c_inner,
                 np.ndarray[np.float_t, ndim=1] c_outer,
                 double fva_tol,
                 int verbosity):
        """ A class to handle the bilevel optimization required by DFVA """

        # Using the old style here since I had issues with cython
        # inheretance...
        GLPKfba.__init__(self, cobra_model, verbosity)

        self.c_inner = c_inner
        self.c_outer = c_outer
        self.fva_tol = fva_tol

        self.lb_orig = np.empty(self.n, dtype=np.float)
        self.ub_orig = np.empty(self.n, dtype=np.float)


    cpdef int solve(self):

        cdef np.float_t[:] temp_fluxes = np.empty(self.n, dtype=np.float)
        cdef int i, flag
        cdef np.float diff

        # Store original bounds prior to the FVA solve
        for i in range(self.n):
            self.lb_orig[i] = self.get_lb_i(i)
            self.ub_orig[i] = self.get_ub_i(i)

        # First solve the inner problem
        for i in range(self.n): self.set_obj_i(i, self.c_inner[i])
        flag = GLPKfba.solve(self)
        if flag == -1: return -1

        # Store resulting fluxes
        for i in range(self.n):
            temp_fluxes[i] = self.get_flux_i(i)

        # Update bounds accordingly
        for i in range(self.n):
            if abs(self.c_inner[i]) >= 1E-6: 
                diff = abs((1. - self.fva_tol) * temp_fluxes[i])
                self.set_bounds_i(
                    i, temp_fluxes[i] - diff, temp_fluxes[i] + diff)

        # Now solve outer problem
        for i in range(self.n): self.set_obj_i(i, self.c_outer[i])
        flag = GLPKfba.solve(self)
        if flag == -1: return -1

        # Now reset bounds to their original values
        for i in range(self.n): 
            self.set_bounds_i(i, self.lb_orig[i], self.ub_orig[i])

        return 0
        



        
