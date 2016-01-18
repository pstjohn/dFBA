import numpy as np
cimport numpy as np

from warnings import warn

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
                 ReactionBoundsFunction bounds_func,
                 collect_fluxes=False,
                 ):


        # Initialize LP solution
        self.fba_prob = glp.GLPKfba(cobra_model)

        # Initialize ODE flags
        self._has_run = 0
        self._has_found_ydot = 0

        # Initialize flux collection if requested
        if collect_fluxes:
            self._collect_fluxes = 1
            self.NV = len(cobra_model.reactions)

        self._y0 = y0 # Store y0 for later re-initialization.

        # Initialize death rate
        self.death_rate = 0.

        self.bounds_func = bounds_func
        assert y0.shape[0] == reaction_indices.shape[0], "y0 and index shape mismatch"
        self.NEQ = y0.shape[0]

        # Initialize external reaction variables
        self.reaction_indices = reaction_indices
        self.ub_e = np.empty(reaction_indices.shape[0])
        self.lb_e = np.empty(reaction_indices.shape[0])
        self.v_e = np.empty(reaction_indices.shape[0])
        for i in range(reaction_indices.shape[0]):
            self.ub_e[i] = cobra_model.reactions[reaction_indices[i]].upper_bound
            self.lb_e[i] = cobra_model.reactions[reaction_indices[i]].lower_bound

        self.y = N_VNew_Serial(self.NEQ)
        nv.arr2nv(self.y, y0)


    def update_bounds_func_parameters(self, np.ndarray[np.float_t, ndim=1]
                                      parameters):
        """ Update the parameters stored in the attached ReactionBoundsFunction
        """

        self.bounds_func.p = parameters


    cdef int _initialize_cvodes(self, realtype t0):
        """ Intialize cvode object """

        # Some problem-specific constants
        cdef int flag
        cdef realtype reltol = 1e-6     # scalar relative tolerance            
        cdef realtype abstol = 1e-8     # scalar relative tolerance            
        cdef long maxnumsteps = 5000

        # Initialize ODE
        flag = CVodeInit(self.cvode_mem, f, t0, self.y)
        flag = CVodeSetUserData(self.cvode_mem, <void*> self)

        flag = CVodeSStolerances(self.cvode_mem, reltol, abstol)
        flag = CVodeSetMaxNumSteps(self.cvode_mem, maxnumsteps)

        flag = CVDense(self.cvode_mem, self.NEQ)

        return flag


    cpdef void reset(self):
        """ Reset the integration to start at the y0 passed during class
        initialization """

        nv.arr2nv(self.y, np.asarray(self._y0))

    cpdef void set_death_rate(self, realtype death_rate):
        """ Set the specific death rate """

        self.death_rate = death_rate

    cpdef int integrate(self, realtype t_start=0., realtype t_end=1., int res=200) except -1:

        cdef int i=0, iout=1

        if self._has_run == 0:
            self._initialize_cvodes(t_start)

        else: 
            CVodeReInit(self.cvode_mem, t_start, self.y)

        # Initialize output arrays
        self._ts = np.empty(res, dtype=np.float)
        self._ys = np.empty([res, self.NEQ], dtype=np.float)
        self._ys_dot = np.empty([res, self.NEQ], dtype=np.float)

        self._ys[:] = NAN
        self._ys_dot[:] = NAN


        if self._collect_fluxes:
            self._vs = np.empty([res, self.NV], dtype=np.float)

        self._ts = np.linspace(t_start, t_end, res)

        # Add initial time point
        for i in range(self.NEQ): self._ys[0,i] = NV_Ith_S(self.y, i)
        self._ts[0] = t_start

        while True:

            flag = CVode(self.cvode_mem, self._ts[iout], self.y, &t_start, CV_NORMAL)

            if flag == CV_SUCCESS:
    
                # Store solution vector
                for i in range(self.NEQ):
                    self._ys[iout, i] = NV_Ith_S(self.y, i)

                # Store optimal flux vector
                if self._collect_fluxes:
                    for i in range(self.NV):
                        self._vs[iout, i] = self.fba_prob.get_flux_i(i)

                iout += 1

            # Break on error flag
            if flag < 0: break

            # Break on finished loop
            elif iout >= res: break

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
        cdef int i = 0, index = 0, flag
        for i in range(self.NEQ):
            index = self.reaction_indices[i]
            self.fba_prob.set_bounds_i(index, self.lb_e[i], self.ub_e[i])

        flag = self.fba_prob.solve()

        for i in range(self.NEQ):
            index = self.reaction_indices[i]
            self.v_e[i] = self.fba_prob.get_flux_i(index)

        return flag

    def _test_fba_solution(self): 
        return self.fba_prob.solve()

    def glp_get_status(self):    return self.fba_prob.get_status()
    def glp_get_prim_stat(self): return self.fba_prob.get_prim_stat()
    def glp_get_dual_stat(self): return self.fba_prob.get_dual_stat()
    def glp_update_from_model(self, cobra_model):
        return self.fba_prob.update_from_model(cobra_model)

    property ys:
        def __get__(self): return np.asarray(self._ys)

    property vs:
        def __get__(self): 
            if self._collect_fluxes:
                return np.asarray(self._vs)
            else: raise RuntimeWarning('Fluxes not collected')

    property ts:
        def __get__(self): return np.asarray(self._ts)

    property ys_dot:
        def __get__(self):

            if self._has_found_ydot == 0:
                self.get_ydot()

            return np.asarray(self._ys_dot)

    cdef int get_ydot(self) except -1:
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



    cpdef int print_final_stats(self) except -1:

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

        return 0


cdef int f(realtype t, N_Vector y_nv, N_Vector ydot_nv, void *user_data):
    
    # Recast variables to more managable forms
    cdef realtype* y = (<N_VectorContent_Serial>y_nv.content).data
    cdef realtype* ydot = (<N_VectorContent_Serial>ydot_nv.content).data
    cdef DFBA_Simulator self = <DFBA_Simulator> user_data

    cdef int flag
    flag = self.find_external_fluxes(y)
    if flag == -1: return -1

    cdef int i = 0
    for i in range(0, self.NEQ):
        ydot[i] = y[0] * self.v_e[i]
        

    # Death rate
    ydot[0] += -self.death_rate * y[0]

    return 0


cdef void check_cvode_flag(int flag):

    if flag == CV_TOO_MUCH_WORK: 
        warn("CVODE: The solver took mxstep internal steps but could not reach tout")
    elif flag == CV_TOO_MUCH_ACC:
        warn("CVODE: The solver could not satisfy the accuracy demanded by the user for some internal step")
    elif flag == CV_ERR_FAILURE:
        warn("CVODE: Error test failures occurred too many times during one internal time step or minimum step size was reached")
    elif flag == CV_CONV_FAILURE:
        warn("CVODE: Convergence test failures occurred too many times during one internal time step or minimum step size was reached")
    elif flag == CV_LINIT_FAIL:
        warn("CVODE: The linear solver's initialization function failed")
    elif flag == CV_LSETUP_FAIL:
        warn("CVODE: The linear solver's setup function failed in an unrecoverable manner")
    elif flag == CV_LSOLVE_FAIL:
        warn("CVODE: The linear solver's solve function failed in an unrecoverable manner")
    elif flag == CV_RHSFUNC_FAIL:
        warn("CVODE: The right-hand side function failed in an unrecoverable manner")
    elif flag == CV_FIRST_RHSFUNC_ERR:
        warn("CVODE: The right-hand side function failed at the first call")
    elif flag == CV_REPTD_RHSFUNC_ERR:
        warn("CVODE: The right-hand side function had repetead recoverable errors")
    elif flag == CV_UNREC_RHSFUNC_ERR:
        warn("CVODE: The right-hand side function had a recoverable error, but no recovery is possible")
    elif flag == CV_RTFUNC_FAIL:
        warn("CVODE: The rootfinding function failed in an unrecoverable manner")
    elif flag == CV_MEM_FAIL:
        warn("CVODE: A memory allocation failed")
    elif flag == CV_MEM_NULL:
        warn("CVODE: The cvode_mem argument was NULL")
    elif flag == CV_ILL_INPUT:
        warn("CVODE: One of the function inputs is illegal")
    elif flag == CV_NO_MALLOC:
        warn("CVODE: The CVODE memory block was not allocated by a call to CVodeMalloc")
    elif flag == CV_BAD_K:
        warn("CVODE: The derivative order k is larger than the order used")
    elif flag == CV_BAD_T:
        warn("CVODE: The time t s outside the last step taken")
    elif flag == CV_BAD_DKY:
        warn("CVODE: The output derivative vector is NULL")
    elif flag == CV_TOO_CLOSE:
        warn("CVODE: The output and initial times are too close to each other")



    # warnings.warn_explicit("a warning", category=UserWarning, filename=__FILE__, lineno=__LINE__)


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
