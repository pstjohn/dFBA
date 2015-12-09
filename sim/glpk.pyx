from libc.stdio cimport printf
import numpy as np
cimport numpy as np


cdef class GLPKfba:

    def __cinit__(self):
        self.lp = glp_create_prob()

    def __dealloc__(self):
        glp_delete_prob(self.lp)
        glp_free_env()

    def __init__(self, 
                 np.ndarray[np.int32_t,   ndim=1, mode="c"] ia, 
                 np.ndarray[np.int32_t,   ndim=1, mode="c"] ja, 
                 np.ndarray[np.float64_t, ndim=1, mode="c"] ar,
                 int nm, 
                 int nr):
        """ Initialize the class with the stochiometric matrix A, and
        initialize the GLP problem """

        # Initialize objective
        glp_set_prob_name(self.lp, "FBA")
        glp_set_obj_dir(self.lp, GLP_MAX) 

        self.m = nm # Metabolites
        self.n = nr # Reactions

        cdef int tot_size = self.m * self.n
        cdef int i, j

        # Allocate space for the GLP solution.
        glp_add_rows(self.lp, self.m)
        glp_add_cols(self.lp, self.n)

        # Add stochiometric constraints
        for i in xrange(self.m): glp_set_row_bnds(self.lp, i+1, GLP_FX, 0.0, 0.0)

        # Load stochiometry into GLP. Numpy array pointers give direct access
        # to the underlying numpy data.
        glp_load_matrix(self.lp, ia.shape[0] - 1, <int*> &ia[0], <int*> &ja[0],
                        <double*> &ar[0])

        # Silence the GLP output
        glp_term_hook(silent_hook, NULL)

    cdef int set_objective(self, np.ndarray[np.float_t, ndim=1] c):

        # Check input argument
        assert c.shape[0] == self.n, "Objective vector length mismatch"

        for j in xrange(self.n):
            glp_set_obj_coef(self.lp, j+1, c[j])

        return 0

    cdef int set_bounds(self, 
                         np.ndarray[np.float_t, ndim=1] lb,
                         np.ndarray[np.float_t, ndim=1] ub):

        # Check input argument
        assert lb.shape[0] == self.n, "Lower-bound vector length mismatch"
        assert ub.shape[0] == self.n, "Upper-bound vector length mismatch"

        cdef int j
        for j in xrange(self.n):
            glp_set_col_bnds(self.lp, j+1, GLP_DB, lb[j], ub[j])

        return 0

    cdef int set_bounds_i(self, int i, float lb, float ub):
        glp_set_col_bnds(self.lp, i+1, GLP_DB, lb, ub)
        return 0

    cdef int solve(self):
        cdef int ret
        ret = glp_simplex(self.lp, NULL)
        if ret != 0: raise RuntimeError("GLP Simplex Failed")

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
