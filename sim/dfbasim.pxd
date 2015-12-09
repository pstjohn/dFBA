import numpy as np
cimport numpy as np

cimport dFBA.sim.glpk as glp
from dFBA.sim.cvode cimport *
from dFBA.sim.ReactionBoundsFunction cimport ReactionBoundsFunction

cdef class DFBA_Simulator:
    cdef:
        void *cvode_mem
        int flag, flagr, iout, NEQ, NOUT, _has_run, _has_found_ydot
        N_Vector y
        glp.GLPKfba fba_prob
        np.int32_t [:] reaction_indices
        double [:] ub_e
        double [:] lb_e
        np.float_t [:] v_e

        np.float_t [:] _y0
        np.float_t [:] _ts
        np.float_t [:, :] _ys
        np.float_t [:, :] _ys_dot

        ReactionBoundsFunction bounds_func

    cpdef int integrate(self, realtype t_out, int res)
    cdef int find_external_fluxes(self, realtype* y)
    cpdef void print_final_stats(self)
    cdef int get_ydot(self)
