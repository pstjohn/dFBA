import numpy as np
cimport numpy as np

cdef class GLPKfba:

    cdef glp_prob *lp
    cdef glp_smcp *opt
    cdef int m, n, verbosity
    cdef int set_objective(self, np.ndarray[np.float_t, ndim=1] c)
    cdef int set_obj_i(self, int j, double obj)
    cdef int set_bounds(self, 
                        np.ndarray[np.float_t, ndim=1] lb,
                        np.ndarray[np.float_t, ndim=1] ub)
    cdef int set_bounds_i(self, int i, float lb, float ub)
    cpdef int solve(self)
    cdef double get_objective(self)
    cpdef np.ndarray[np.float_t, ndim=1] get_fluxes(self)
    cdef float get_flux_i(self, int i)
    cdef double get_lb_i(self, int i)
    cdef double get_ub_i(self, int i)
    

    cdef int get_status(self)
    cdef int get_prim_stat(self)
    cdef int get_dual_stat(self)
    

cdef extern from "glpk.h":
    ctypedef struct glp_prob
    ctypedef struct glp_smcp:
        int msg_lev
        int meth
        int presolve
        double tol_bnd
        double tol_dj
        double obj_ul

    void glp_set_prob_name(glp_prob *P, const char *name)
    void glp_set_obj_dir(glp_prob *P, int dir)
    int glp_add_rows(glp_prob *P, int nrs)
    void glp_set_row_name(glp_prob *P, int i, const char *name)
    void glp_set_row_bnds(glp_prob *P, int i, int type, double lb, double ub)
    int glp_add_cols(glp_prob *P, int ncs)
    void glp_set_col_name(glp_prob *P, int j, const char *name)
    void glp_set_col_bnds(glp_prob *P, int j, int type, double lb, double ub)
    void glp_set_obj_coef(glp_prob *P, int j, double coef)
    void glp_load_matrix(glp_prob *P, int ne, const int ia[], const int ja[],
                         const double ar[])
    int glp_simplex(glp_prob *P, const glp_smcp *parm)
    double glp_get_obj_val(glp_prob *P)
    double glp_get_col_prim(glp_prob *P, int j)
    double glp_get_col_dual(glp_prob *P, int j)
    glp_prob *glp_create_prob()
    void glp_delete_prob(glp_prob *P)
    int glp_free_env()
    void glp_term_hook(int (*func)(void *info, const char *s), void *info)
    int glp_get_status(glp_prob *P)
    int glp_get_prim_stat(glp_prob *P)
    int glp_get_dual_stat(glp_prob *P)

    double glp_get_col_lb(glp_prob *P, int j);
    double glp_get_col_ub(glp_prob *P, int j);



    enum: GLP_MAX, GLP_UP, GLP_LO, GLP_MSG_OFF, GLP_FX, GLP_DB

    enum: GLP_EBADB, GLP_ESING, GLP_ECOND, GLP_EBOUND, GLP_EFAIL, GLP_EOBJLL
    enum: GLP_EOBJUL, GLP_EITLIM, GLP_ETMLIM, GLP_ENOPFS, GLP_ENODFS

    enum: GLP_UNDEF, GLP_FEAS, GLP_INFEAS, GLP_NOFEAS, GLP_OPT, GLP_UNBND

