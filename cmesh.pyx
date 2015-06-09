
cimport cython

cdef extern from "mesh.h":
    ctypedef struct Mesh:
        int nn, ne, nl
        double *x
        double *y
        int *bnd
        int *ele
        int *neigh
        int *edge
        int *bnd_edge

    Mesh ReadMesh(char *filename)
    void FreeMesh(Mesh *m)

cdef class CMesh:
    cdef Mesh m

    def __init__(self, filename):
        cdef bytes py_bytes = filename.encode()
        cdef char *c_filename = py_bytes
        m = ReadMesh(c_filename)