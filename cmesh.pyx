cimport cython
cimport numpy
import numpy

from libc.stdlib cimport malloc

cdef extern from "mesh.h":
    ctypedef struct Mesh:
        int nn, ne, nl
        double* x
        double* y
        int *bnd
        int *ele
        int *neigh
        int *edge
        int *bnd_edge

    Mesh ReadMesh(char *filename)
    void FreeMesh(Mesh *m)

cdef class CMesh:
    cdef Mesh* m

    def __init__(self, filename):
        cdef bytes py_bytes = filename.encode()
        cdef char *c_filename = py_bytes
        self.m = <Mesh*>malloc(sizeof(Mesh))
        self.m[0] = ReadMesh(c_filename)

    def __dealloc__(self):
        FreeMesh(self.m)

    @property
    def x(self):
        cdef int n = self.m.nn
        cdef double[:] x = <double[:n]>self.m.x
        return numpy.array(x, dtype=numpy.double)
