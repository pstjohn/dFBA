import numpy as np
cimport numpy as np


# Functions to convert to/from NVector & DlsMat objects to Numpy Arrays.
# These routines are pass-by-reference, 

cdef void nv2mem_view(N_Vector v, np.float_t [:] p):
    
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    assert n == p.shape[0], "vector length mismatch"

    cdef long int i
    for i in range(n):
        p[i] = v_data[i]

cdef void mem_view2nv(N_Vector v, np.float_t [:] p):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    assert n == p.shape[0], "vector length mismatch"

    cdef long int i
    for i in range(n):
        v_data[i] = p[i]

cdef void nv2arr(N_Vector v, np.ndarray[np.float_t, ndim=1] p):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    assert n == p.shape[0], "vector length mismatch"

    cdef long int i
    for i in range(n):
        p[i] = v_data[i]

cdef void arr2nv(N_Vector v, np.ndarray[np.float_t, ndim=1] p):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    assert n == p.shape[0], "vector length mismatch"

    cdef long int i
    for i in range(n):
        v_data[i] = p[i]

cdef void dls2np(DlsMat A, np.ndarray[np.float_t, ndim=2] P):
    cdef long int n = A.N
    cdef long int m = A.M
    
    assert m == P.shape[0], "row shape mismatch"
    assert n == P.shape[1], "col shape mismatch"

    cdef long int i, j
    cdef realtype* col_j

    # Iteration scheme as recommended by sundial's header.
    for j in range(n):
        col_j = A.cols[j]
        for i in range(m):
            P[i,j] = col_j[i]


cdef void np2dls(DlsMat A, np.ndarray[np.float_t, ndim=2] P):
    cdef long int n = A.N
    cdef long int m = A.M
    
    assert m == P.shape[0], "row shape mismatch"
    assert n == P.shape[1], "col shape mismatch"

    cdef long int i, j
    cdef realtype* col_j

    # Iteration scheme as recommended by sundial's header.
    for j in range(n):
        col_j = A.cols[j]
        for i in range(m):
            col_j[i] = P[i,j]

