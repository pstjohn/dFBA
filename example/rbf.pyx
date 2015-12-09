from dFBA.sim.dfbasim cimport ReactionBoundsFunction, realtype

cdef class GlucoseUptake(ReactionBoundsFunction):
    cdef int evaluate(self, realtype* y, 
                      double [:] lb_e, 
                      double [:] ub_e): 

        # Glucose MM uptake.
        lb_e[1] = -10 * y[1] / (0.5 + y[1])

        return 0
