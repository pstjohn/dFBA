import numpy as np
cimport numpy as np

from dFBA.sim.dfbasim cimport DFBA_Simulator 
from dFBA.sim.glpk_fva cimport GLPKfva
from dFBA.sim.ReactionBoundsFunction cimport ReactionBoundsFunction


cdef class DFVA_Simulator_base(DFBA_Simulator):
    cdef:
        np.float_t[:] c_inner
        np.float_t[:] c_outer
        double fva_tol

