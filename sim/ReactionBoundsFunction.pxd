from dFBA.sim.cvode cimport realtype

cdef class ReactionBoundsFunction:
    cdef double [:] p # Parameter array
    cdef int evaluate(self, realtype* y, double [:] lb_e, double [:] ub_e)


