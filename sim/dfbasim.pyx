import numpy as np
cimport numpy as np

cimport dFBA.sim.nvector as nv
from dFBA.sim.ReactionBoundsFunction cimport ReactionBoundsFunction

cdef class DFBA_Simulator:
    # cdef:
    #     void *cvode_mem
    #     int flag, flagr, iout, NEQ, NOUT, _has_run
    #     N_Vector y
    #     glp.GLPKfba fba_prob
    #     np.int32_t [:] reaction_indices
    #     np.float_t [:] ub_e
    #     np.float_t [:] lb_e
    #     np.float_t [:] v_e

    #     np.float_t [:] _y0
    #     np.float_t [:] _ts
    #     np.float_t [:, :] _ys

    #     ReactionBoundsFunction bounds_func

    def __cinit__(self):
        self.cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON)


    def __dealloc__(self):
        if self.y: N_VDestroy_Serial(self.y)
        CVodeFree(&self.cvode_mem)


    def __init__(self, cobra_model, 
                 np.ndarray[np.int32_t, ndim=1] reaction_indices,
                 np.ndarray[np.float_t, ndim=1] y0, 
                 ReactionBoundsFunction bounds_func
                 ):


        ## Initialize Linear Programming Object.
        # Should I move this to the GLPK class?

        # Parse cobra_model
        cdef int nm, nr, i, j, tot_size
        cdef float stoich
        nm = len(cobra_model.metabolites)
        nr = len(cobra_model.reactions)

        # Initialize objective, bounds, and stoich containters
        cdef np.ndarray[np.float_t, ndim=1] c  = np.empty(nr, dtype=np.float)
        cdef np.ndarray[np.float_t, ndim=1] lb = np.empty(nr, dtype=np.float)
        cdef np.ndarray[np.float_t, ndim=1] ub = np.empty(nr, dtype=np.float)

        # glpk wants the indexing to start at 1 for some reason.
        ia_list = [0.]
        ja_list = [0.]
        ar_list = [0.]

        # Fill containers
        for j in xrange(nr):
            reaction = cobra_model.reactions[j]
            c[j] = reaction.objective_coefficient
            ub[j] = reaction.upper_bound
            lb[j] = reaction.lower_bound
            for metabolite, stoich in reaction.metabolites.iteritems():
                i = cobra_model.metabolites.index(metabolite.id)
                ia_list += [i+1]
                ja_list += [j+1]
                ar_list += [stoich]

        # Allocate+fill sparse numpy arrays for stoich
        tot_size = len(ia_list)

        cdef np.ndarray[np.int32_t,   ndim=1, mode="c"] ia = \
            np.zeros(tot_size, dtype=np.int32)
        cdef np.ndarray[np.int32_t,   ndim=1, mode="c"] ja = \
            np.zeros(tot_size, dtype=np.int32)
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] ar = \
            np.zeros(tot_size, dtype=np.float64)

        for i in xrange(tot_size):
            ia[i] = ia_list[i]
            ja[i] = ja_list[i]
            ar[i] = ar_list[i]

        # Initialize LP solution
        self.fba_prob = glp.GLPKfba(ia, ja, ar, nm, nr)
        self.fba_prob.set_bounds(lb, ub)
        self.fba_prob.set_objective(c)
        

        # Initialize ODE flags
        self._has_run = 0
        self._has_found_ydot = 0

        self._y0 = y0 # Store y0 for later re-initialization.

        self.bounds_func = bounds_func
        assert y0.shape[0] == reaction_indices.shape[0], "y0 and index shape mismatch"
        self.NEQ = y0.shape[0]

        # Initialize external reaction variables
        self.reaction_indices = reaction_indices
        self.ub_e = np.empty(reaction_indices.shape[0])
        self.lb_e = np.empty(reaction_indices.shape[0])
        self.v_e = np.empty(reaction_indices.shape[0])
        for i in range(reaction_indices.shape[0]):
            self.ub_e[i] = ub[reaction_indices[i]]
            self.lb_e[i] = lb[reaction_indices[i]]

        self.y = N_VNew_Serial(self.NEQ)
        nv.arr2nv(self.y, y0)

        # Some problem-specific constants
        cdef realtype t0     = 0.0      # initial time           
        cdef realtype reltol = 1e-6     # scalar relative tolerance            
        cdef realtype abstol = 1e-8     # scalar relative tolerance            
        cdef long maxnumsteps = 5000

        # Initialize ODE
        flag = CVodeInit(self.cvode_mem, f, t0, self.y)
        flag = CVodeSetUserData(self.cvode_mem, <void*> self)

        flag = CVodeSStolerances(self.cvode_mem, reltol, abstol)
        flag = CVodeSetMaxNumSteps(self.cvode_mem, maxnumsteps)


        flag = CVDense(self.cvode_mem, self.NEQ)
        # flag = CVDlsSetDenseJacFn(self.cvode_mem, Jac)

    cpdef int integrate(self, realtype t_out, int res):

        cdef realtype t = 0.
        cdef realtype tnext = 0.
        cdef int i=0, iout=1

        if self._has_run == 1:
            nv.arr2nv(self.y, np.asarray(self._y0))
            CVodeReInit(self.cvode_mem, t, self.y)

        # Initialize output arrays
        self._ts = np.empty(res, dtype=np.float)
        self._ys = np.empty([res, self.NEQ], dtype=np.float)
        self._ys_dot = np.empty([res, self.NEQ], dtype=np.float)

        self._ts = np.linspace(t, t_out, res)

        # Add initial time point
        for i in range(self.NEQ): self._ys[0,i] = NV_Ith_S(self.y, i)
        self._ts[0] = t

        while True:

            flag = CVode(self.cvode_mem, self._ts[iout], self.y, &t, CV_NORMAL)

            if flag == CV_SUCCESS:
                for i in range(self.NEQ):
                    self._ys[iout, i] = NV_Ith_S(self.y, i)
                iout += 1

            else: raise RuntimeError("CVode Exception")

            if iout >= res: break

        self._has_run = 1
        self._has_found_ydot = 0
        return 0



    cdef int find_external_fluxes(self, realtype* y):
        """ Calculate the optimal external exchange fluxes for the given y, and
        return the fluxes in self.v_e. Makes use of the *hopefully* external
        external_reaction_bounds """

        # Use external enzyme chemistry to find new upper and lower bounds
        self.bounds_func.evaluate(y, self.lb_e, self.ub_e)

        # Propogate external bounds to linear program
        cdef int i = 0, index = 0
        for i in range(self.NEQ):
            index = self.reaction_indices[i]
            self.fba_prob.set_bounds_i(index, self.lb_e[i], self.ub_e[i])

        self.fba_prob.solve()

        for i in range(self.NEQ):
            index = self.reaction_indices[i]
            self.v_e[i] = self.fba_prob.get_flux_i(index)

        return 0

    property ys:
        def __get__(self): return np.asarray(self._ys)

    property ts:
        def __get__(self): return np.asarray(self._ts)

    property ys_dot:
        def __get__(self):

            if self._has_found_ydot == 0:
                self.get_ydot()

            return np.asarray(self._ys_dot)

    cdef int get_ydot(self):
        """ Get the derivative at each time step """
        assert self._has_run == 1, "Integration has not yet run"

        cdef int i
        cdef N_Vector yi
        cdef N_Vector yi_dot
        yi = N_VNew_Serial(self.NEQ)
        yi_dot = N_VNew_Serial(self.NEQ)

        try:

            for i in range(self._ts.shape[0]):

                # Call function derivative for the current yi
                nv.mem_view2nv(yi, self._ys[i,:])
                f(<realtype> self._ts[i], yi, yi_dot, <void*> self)
                nv.nv2mem_view(yi_dot, self._ys_dot[i,:])

        finally:
            # Clean up N_Vectors
            N_VDestroy_Serial(yi)
            N_VDestroy_Serial(yi_dot)

        return 0 



    cpdef void print_final_stats(self):

        cdef long nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf

        CVodeGetNumSteps(self.cvode_mem, &nst)
        CVodeGetNumRhsEvals(self.cvode_mem, &nfe)
        CVodeGetNumLinSolvSetups(self.cvode_mem, &nsetups)
        CVodeGetNumErrTestFails(self.cvode_mem, &netf)
        CVodeGetNumNonlinSolvIters(self.cvode_mem, &nni)
        CVodeGetNumNonlinSolvConvFails(self.cvode_mem, &ncfn)

        CVDlsGetNumJacEvals(self.cvode_mem, &nje)
        CVDlsGetNumRhsEvals(self.cvode_mem, &nfeLS)

        print "\nFinal Statistics:"
        print "nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld" % (
            nst, nfe, nsetups, nfeLS, nje)
        print "nni = %-6ld ncfn = %-6ld netf = %-6ld\n" %  (nni, ncfn, netf)


cdef int f(realtype t, N_Vector y_nv, N_Vector ydot_nv, void *user_data):
    
    # Recast variables to more managable forms
    cdef realtype* y = (<N_VectorContent_Serial>y_nv.content).data
    cdef realtype* ydot = (<N_VectorContent_Serial>ydot_nv.content).data
    cdef DFBA_Simulator self = <DFBA_Simulator> user_data

    self.find_external_fluxes(y)

    # Biomass growth
    ydot[0] = y[0] * self.v_e[0]

    cdef int i = 0
    for i in range(1, self.NEQ):
        ydot[i] = y[0] * self.v_e[i]
        

    return 0


# cdef class ReactionBoundsFunction:
#     cdef int evaluate(self, realtype* y, 
#                       double [:] lb_e, 
#                       double [:] ub_e): 
#         """ A user-defined function which updates the lower and upper bound of the
#         metabolic exchange fluxes.
        
#         The current external concentrations, y, are passed as inputs. lb_e and ub_e
#         are passed with their current (default) values, so only active bounds need
#         to be changed.
#         """

#         return 0
    




# cdef int Jac(long int N, realtype t, N_Vector y_nv, N_Vector fy_nv, DlsMat J_dls, 
#              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):

#     cdef realtype* y = (<N_VectorContent_Serial>y_nv.content).data
#     cdef realtype** J = J_dls.cols

#     J[0][0] = 1.
#     J[1][1] = 1.
#     J[2][2] = 1E-1 * y[1]
#     J[2][1] = 1E-1 * y[2]

#     return 0




# def main(cobra_model, np.ndarray[np.int32_t, ndim=1] reaction_indices,
#          np.ndarray[np.float_t, ndim=1] y0):

#     cdef DFBA_Simulator dfbasim = DFBA_Simulator(cobra_model, reaction_indices,
#                                                  y0, GlucoseUptake())
#     print 'Initialization successful'



#     dfbasim.fba_prob.solve()
#     cdef float opt
#     opt = dfbasim.fba_prob.get_objective()
#     print "Objective: ",
#     print opt

#     dfbasim.integrate(10.0, 20)
#     print np.asarray(dfbasim._ys)

#     dfbasim.print_final_stats()
#     return 0.
