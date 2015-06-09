
cimport cython
cimport numpy as np
import numpy as np

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
    cdef Mesh *m

    def __init__(self, filename):
        cdef bytes py_bytes = filename.encode()
        cdef char *c_filename = py_bytes
        cdef Mesh m = ReadMesh(c_filename)
        self.m = &m

    def __dealloc__(self):
        FreeMesh(self.m)