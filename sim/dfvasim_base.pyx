import numpy as np

cdef class DFVA_Simulator_base(DFBA_Simulator):

    def __init__(self, cobra_model, 
                 np.ndarray[np.int32_t, ndim=1] reaction_indices,
                 np.ndarray[np.float_t, ndim=1] y0, 
                 ReactionBoundsFunction bounds_func,
                 np.ndarray[np.float_t, ndim=1] c_inner,
                 np.ndarray[np.float_t, ndim=1] c_outer,
                 double fva_tol=.99,
                 ):

        self.c_inner = c_inner
        self.c_outer = c_outer
        self.fva_tol = fva_tol

        DFBA_Simulator.__init__(self, cobra_model, reaction_indices, y0,
                                bounds_func)
    

    def _initialize_lp(self, cobra_model):
        """ Initailize the fva lp """

        cdef np.ndarray[np.float_t, ndim=1] c_inner = np.asarray(self.c_inner)
        cdef np.ndarray[np.float_t, ndim=1] c_outer = np.asarray(self.c_outer)
        self.fba_prob = GLPKfva(cobra_model, c_inner, c_outer, self.fva_tol)

    


