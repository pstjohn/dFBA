cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector
    N_Vector N_VNew_Serial(long vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)
    int NV_LENGTH_S(N_Vector v)
    realtype NV_Ith_S(N_Vector v, int i)

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat: 
        int type 
        long M 
        long N 
        long ldim 
        long mu 
        long ml 
        long s_mu 
        realtype *data 
        long ldata 
        realtype **cols

    ctypedef _DlsMat* DlsMat

    void PrintMat(DlsMat A)
    DlsMat NewDenseMat(long M, long N)
    void DestroyMat(DlsMat A)

cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long length
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial

cdef extern from "cvode/cvode.h":
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data)
    ctypedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data)
    ctypedef int (*CVDlsDenseJacFn)(long N, realtype t, 
             N_Vector y, N_Vector fy,  
             DlsMat Jac, void *user_data, 
             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    ctypedef void (*CVErrHandlerFn)(int error_code, const char *module, const
                                    char *function, char *msg, void *eh_data)

    # Core solver functions
    void* CVodeCreate(int lmm, int iter)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret,
              int itask)
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
    void CVodeFree(void **cvode_mem)

    # Solver options
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetUserData(void *cvode_mem, void *user_data)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac)
    int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void *eh_data)
    int CVodeSetMinStep(void *cvode_mem, realtype hmin)

    # Solver Stats Functions
    int CVodeGetNumSteps(void *cvode_mem, long *nsteps)
    int CVodeGetNumRhsEvals(void *cvode_mem, long *nfevals)
    int CVodeGetNumLinSolvSetups(void *cvode_mem, long *nlinsetups) 
    int CVodeGetNumErrTestFails(void *cvode_mem, long *netfails)
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long *nniters)
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long *nncfails)
   

    enum: CV_ADAMS, CV_BDF, CV_NEWTON

    enum: CV_NORMAL, CV_SUCCESS

    enum: CV_TOO_MUCH_ACC, CV_ERR_FAILURE, CV_CONV_FAILURE, CV_LINIT_FAIL
    enum: CV_TOO_MUCH_WORK, CV_LSETUP_FAIL, CV_LSOLVE_FAIL, CV_RHSFUNC_FAIL
    enum: CV_FIRST_RHSFUNC_ERR, CV_REPTD_RHSFUNC_ERR, CV_UNREC_RHSFUNC_ERR
    enum: CV_RTFUNC_FAIL, CV_MEM_FAIL, CV_MEM_NULL, CV_ILL_INPUT
    enum: CV_NO_MALLOC, CV_BAD_K, CV_BAD_T, CV_BAD_DKY, CV_TOO_CLOSE

    
    
cdef extern from "cvode/cvode_direct.h":
    int CVDlsGetNumJacEvals(void *cvode_mem, long *njevals)
    int CVDlsGetNumRhsEvals(void *cvode_mem, long *nfevalsLS)

cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, long N)

