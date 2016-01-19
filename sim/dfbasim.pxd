import numpy as np
cimport numpy as np

cimport dFBA.sim.glpk as glp
from dFBA.sim.cvode cimport *
from dFBA.sim.ReactionBoundsFunction cimport ReactionBoundsFunction
cimport dFBA.sim.nvector as nv

from warnings import warn

cdef class DFBA_Simulator:
    cdef:
        void *cvode_mem
        int flag, flagr, iout, NEQ, NOUT, _has_run, _has_found_ydot,
        int _collect_fluxes, NV
        N_Vector y
        glp.GLPKfba fba_prob
        np.int32_t [:] reaction_indices
        double [:] ub_e
        double [:] lb_e
        np.float_t [:] v_e
        realtype death_rate

        np.float_t [:] _y0
        np.float_t [:] _ts
        np.float_t [:, :] _ys
        np.float_t [:, :] _ys_dot
        np.float_t [:, :] _vs

        ReactionBoundsFunction bounds_func

    cpdef int integrate(self, realtype t_start=*, realtype t_end=*, int res=*) except -1
    cdef int find_external_fluxes(self, realtype* y)
    cpdef int print_final_stats(self) except -1
    cdef int get_ydot(self) except -1
    cdef int get_vs(self) except -1
    cpdef void reset(self)
    cdef int _initialize_cvodes(self, realtype t0)
    cpdef void set_death_rate(self, realtype death_rate)
    
    
    

cdef extern from 'math.h':
    float NAN
