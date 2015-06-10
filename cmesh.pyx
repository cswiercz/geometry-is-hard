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
        cdef int nn = self.m.nn
        cdef double[:] x = <double[:nn]>self.m.x
        return numpy.array(x, dtype=numpy.double)

    @property
    def y(self):
        cdef int nn = self.m.nn
        cdef double[:] y = <double[:nn]>self.m.y
        return numpy.array(y, dtype=numpy.double)

    @property
    def triangles(self):
        cdef int ne = self.m.ne
        cdef int[:,:] ele = <int[:ne,:3]>self.m.ele
        return numpy.array(ele, dtype=numpy.int)

    @property
    def edges(self):
        cdef int nl = self.m.nl
        cdef int[:,:] edge = <int[:nl,:2]>self.m.edge
        return numpy.array(edge, dtype=numpy.int)

    @property
    def neighbors(self):
        cdef int ne = self.m.ne
        cdef int[:,:] neigh = <int[:ne,:3]>self.m.neigh
        return numpy.array(neigh, dtype=numpy.int)
