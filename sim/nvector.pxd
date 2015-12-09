import numpy as np
cimport numpy as np

cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat: 
        int type 
        long int M 
        long int N 
        long int ldim 
        long int mu 
        long int ml 
        long int s_mu 
        realtype *data 
        long int ldata 
        realtype **cols

    ctypedef _DlsMat* DlsMat

    void PrintMat(DlsMat A);
    DlsMat NewDenseMat(long int M, long int N);
    void DestroyMat(DlsMat A);

cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial

cdef void nv2arr(N_Vector v, np.ndarray[np.float_t, ndim=1] p)
cdef void nv2mem_view(N_Vector v, np.float_t [:] p)
cdef void mem_view2nv(N_Vector v, np.float_t [:] p)
cdef void arr2nv(N_Vector v, np.ndarray[np.float_t, ndim=1] p)
cdef void dls2np(DlsMat A, np.ndarray[np.float_t, ndim=2] P)
cdef void np2dls(DlsMat A, np.ndarray[np.float_t, ndim=2] P)

