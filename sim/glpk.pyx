from libc.stdio cimport printf
import numpy as np
cimport numpy as np


cdef class GLPKfba:

    def __cinit__(self):
        self.lp = glp_create_prob()

    def __dealloc__(self):
        glp_delete_prob(self.lp)
        # glp_free_env()

    def __init__(self, cobra_model):
        """ Initialize the class with the stochiometric matrix A, and
        initialize the GLP problem """

        # Parse cobra_model
        cdef int nm, nr, i, j
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

        # Initialize objective
        glp_set_prob_name(self.lp, "FBA")
        glp_set_obj_dir(self.lp, GLP_MAX) 

        self.m = nm # Metabolites
        self.n = nr # Reactions

        # Allocate space for the GLP solution.
        glp_add_rows(self.lp, self.m)
        glp_add_cols(self.lp, self.n)

        # Initialize bounds and objective
        self.set_objective(c)
        self.set_bounds(lb, ub)

        # Add stochiometric constraints
        for i in xrange(self.m): glp_set_row_bnds(self.lp, i+1, GLP_FX, 0.0, 0.0)

        # Load stochiometry into GLP. Numpy array pointers give direct access
        # to the underlying numpy data.
        glp_load_matrix(self.lp, ia.shape[0] - 1, <int*> &ia[0], <int*> &ja[0],
                        <double*> &ar[0])

        # Silence the GLP output
        glp_term_hook(silent_hook, NULL)

    cdef int set_objective(self, np.ndarray[np.float_t, ndim=1] c) except -1:

        # Check input argument
        assert c.shape[0] == self.n, "Objective vector length mismatch"

        for j in xrange(self.n):
            glp_set_obj_coef(self.lp, j+1, c[j])

        return 0

    cdef int set_bounds(self, 
                        np.ndarray[np.float_t, ndim=1] lb,
                        np.ndarray[np.float_t, ndim=1] ub) except -1:

        # Check input argument
        assert lb.shape[0] == self.n, "Lower-bound vector length mismatch"
        assert ub.shape[0] == self.n, "Upper-bound vector length mismatch"

        cdef int j
        for j in xrange(self.n):
            if lb[j] == ub[j]:
                glp_set_col_bnds(self.lp, j+1, GLP_FX, lb[j], ub[j])
            elif lb[j] < ub[j]:
                glp_set_col_bnds(self.lp, j+1, GLP_DB, lb[j], ub[j])
            else:
                raise RuntimeError("GLP: Reaction bounds infeasible")

        return 0

    cdef int set_bounds_i(self, int i, float lb, float ub) except -1:
        glp_set_col_bnds(self.lp, i+1, GLP_DB, lb, ub)
        return 0

    cdef int solve(self) except -1:

        cdef int flag
        flag = glp_simplex(self.lp, NULL)
        check_simplex_flag(flag)

        return 0

    cdef double get_objective(self): return glp_get_obj_val(self.lp)

    cdef np.float_t [:] get_fluxes(self):

        cdef np.ndarray[np.float_t, ndim=1] col_prim = np.empty(self.n, dtype=np.float)
        cdef int j

        for j in xrange(self.n): col_prim[j] = glp_get_col_prim(self.lp, j+1)
        return col_prim

    cdef float get_flux_i(self, int i): return glp_get_col_prim(self.lp, i+1)


            

            
            



def main(np.ndarray[np.float_t, ndim=2] A,  # Reduced stochiometric matrix
         np.ndarray[np.float_t, ndim=1] c,  # Objective vector
         np.ndarray[np.float_t, ndim=1] ub, # Upper bound vector
         np.ndarray[np.float_t, ndim=1] lb, # Lower bound vector
         int nr,                            # Number of reactions
         int nm,                            # Number of metabolites
         ):


    cdef GLPKfba glp_prob = GLPKfba(A)
    glp_prob.set_bounds(lb, ub)
    glp_prob.set_objective(c)
    glp_prob.solve()
    print "Objective: ", glp_prob.get_objective()
    print "Fluxes: ", glp_prob.get_fluxes_primary()


    return 0
    

cdef int silent_hook(void *info, const char *s):
    """function to print nothing but trick GLPK into thinking we did"""
    return 1



cdef int check_simplex_flag(int flag) except -1:
    if flag == 0: return 1

    elif flag == GLP_EBADB:
        raise RuntimeError('GLP: Unable to start the search, because the initial basis speci- fied in the problem object is invalid—the number of basic (auxiliary and structural) variables is not the same as the number of rows in the problem object.')
    elif flag == GLP_ESING:
        raise RuntimeError('GLP: Unable to start the search, because the basis matrix corresponding to the initial basis is singular within the working precision.')
    elif flag == GLP_ECOND:
        raise RuntimeError('GLP: Unable to start the search, because the basis matrix corresponding to the initial basis is ill-conditioned, i.e. its condition number is too large.')
    elif flag == GLP_EBOUND:
        raise RuntimeError('GLP: Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds.')
    elif flag == GLP_EFAIL:
        raise RuntimeWarning('GLP: The search was prematurely terminated due to the solver failure.')
    elif flag == GLP_EOBJLL:
        raise RuntimeWarning('GLP: The search was prematurely terminated, because the objective function being maximized has reached its lower limit and continues decreasing (the dual simplex only).')
    elif flag == GLP_EOBJUL:
        raise RuntimeWarning('GLP: The search was prematurely terminated, because the objective function being minimized has reached its upper limit and continues increasing (the dual simplex only).')
    elif flag == GLP_EITLIM:
        raise RuntimeWarning('GLP: The search was prematurely terminated, because the simplex iteration limit has been exceeded.')
    elif flag == GLP_ETMLIM:
        raise RuntimeWarning('GLP: The search was prematurely terminated, because the time limit has been exceeded.')
    elif flag == GLP_ENOPFS:
        raise RuntimeWarning('GLP: The LP problem instance has no primal feasible solution (only elif the LP presolver is used).')
    elif flag == GLP_ENODFS:
        raise RuntimeWarning('GLP: The LP problem instance has no dual feasible solution (only elif the LP presolver is used).')


        # # /* declare variables */

        # cdef int tot_size = A.shape[0]*A.shape[1]
        # # cdef int tot_size = 5

    # cdef np.ndarray[np.int32_t, ndim=1, mode="c"] ia = np.empty(tot_size, dtype=np.int32)
    # cdef np.ndarray[np.int32_t, ndim=1, mode="c"] ja = np.empty(tot_size, dtype=np.int32)
    # cdef np.ndarray[np.float64_t, ndim=1, mode="c"] ar = np.empty(tot_size, dtype=np.float64)

    # cdef double z, x1, x2
    # cdef int i, j, n_elements=0
    # # /* create problem */
    # lp = glp_create_prob()
    # glp_set_prob_name(lp, "short")
    # glp_set_obj_dir(lp, GLP_MAX)
    # # /* fill problem */
    
    # # Add row variables (mass balances)
    # # For now, (S*v = 0), all rows fixed at 0
    # glp_add_rows(lp, A.shape[0])
    # for i in xrange(A.shape[0]):
        # glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0)

    # # Add column variables (reaction rates)
    # # These are defined by the ub and lb vectors
    # glp_add_cols(lp, A.shape[1])
    # for j in xrange(A.shape[1]):
        # glp_set_col_bnds(lp, j+1, GLP_DB, lb[j], ub[j])

    # # Set the objective function
    # for j in xrange(c.shape[0]):
        # glp_set_obj_coef(lp, j+1, c[j])

    # # Get a sparse matrix representation of the A matrix.
    # cdef double tol = 1E-8
    # for i in xrange(A.shape[0]):
        # for j in xrange(A.shape[1]):
        #     if abs(A[i,j]) > tol:
        #         n_elements += 1

        #         ia[n_elements] = i+1
        #         ja[n_elements] = j+1
        #         ar[n_elements] = A[i,j]

    # # for i in xrange(n_elements):
    # #     print ia[i], ja[i], ar[i]
    
    # glp_load_matrix(lp, n_elements, <int*> &ia[0], <int*> &ja[0], <double*> &ar[0])
    # # /* solve problem */
    # glp_term_hook(silent_hook, NULL)
    # glp_simplex(lp, NULL)
    # # /* recover and display results */
    # z = glp_get_obj_val(lp)
    # x1 = glp_get_col_prim(lp, 1)
    # x2 = glp_get_col_prim(lp, 2)
    # printf("z = %g x1 = %g x2 = %g\n", z, x1, x2)
    # # /* housekeeping */
    # glp_delete_prob(lp)
    # glp_free_env()
