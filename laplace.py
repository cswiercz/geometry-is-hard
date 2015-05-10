
import numpy as np
from matplotlib.tri import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


# ------------------------------
def read_triangulation(filename):
    """
    Read in the Triangle mesh stored in <filename>.node, <filename>.ele
    """

    # Read in the vertices of the triangulation
    fid = open(filename + ".node", "r")
    nn = int(fid.readline().split()[0])
    x = np.zeros(nn, dtype = np.float64)
    y = np.zeros(nn, dtype = np.float64)
    bnd = np.zeros(nn, dtype = np.int32)
    for i in range(nn):
        line = fid.readline().split()
        x[i], y[i] = float(line[1]), float(line[2])
        bnd[i] = int(line[3])

    fid.close()

    # Read in the elements of the triangulation
    fid = open(filename + ".ele", "r")
    ne = int(fid.readline().split()[0])
    ele = np.zeros((ne, 3), dtype = np.int32)
    for i in range(ne):
        ele[i,:] = map(int, fid.readline().split()[1:])
    ele -= 1
    fid.close()

    return Triangulation(x, y, ele), bnd


# ----------------------
def mesh_to_matrix(mesh):
    """
    Build a sparse matrix representing the connectivity structure of the
    input mesh, for filling

    Arguments:
    =========
    mesh: a matplotlib.tri.Triangulation object

    Returns:
    =======
    indices, indptr: the internal structure of the CSR matrix object
    """

    nn = len(mesh.x)
    nl, _ = np.shape(mesh.edges)

    indptr = np.ones(nn + 1, dtype = np.int32)
    indices = -np.ones(nn + 2*nl, dtype = np.int32)

    for n in range(nl):
        i, j = mesh.edges[n,:]
        indptr[i+1] += 1
        indptr[j+1] += 1

    indptr[0] = 0
    for n in range(nn):
        indptr[n+1] += indptr[n]

    for i in range(nn):
        indices[indptr[i]] = i

    for n in range(nl):
        i, j = mesh.edges[n,:]
        k = np.argmin(indices[indptr[i] : indptr[i+1]])
        indices[indptr[i] + k] = j

        k = np.argmin(indices[indptr[j] : indptr[j+1]])
        indices[indptr[j] + k] = i

    return indices, indptr


# -------------------
def fe_matrices(mesh):
    """
    Return the stiffness/mass matrices for the Poisson problem on a given
    mesh
    """
    x = mesh.x
    y = mesh.y

    nn = len(x)
    nl, _ = np.shape(mesh.edges)

    indices, indptr = mesh_to_matrix(mesh)

    A = csr_matrix((np.zeros(nn+2*nl, dtype=np.float64), indices, indptr),
                   shape = (nn, nn))
    B = csr_matrix((np.zeros(nn+2*nl, dtype=np.float64), indices, indptr),
                   shape = (nn, nn))

    AE = np.zeros((3, 3))
    BE = np.zeros((3, 3))
    V  = np.zeros((3, 2))

    ne, _ = np.shape(mesh.triangles)
    for n in range(ne):
        ele = mesh.triangles[n,:]

        for i in range(3):
            j = ele[(i+1) % 3]
            k = ele[(i+2) % 3]

            V[i, 0] = y[j] - y[k]
            V[i, 1] = x[k] - x[j]

        area = 0.5 * abs(V[0, 0]*V[1, 1] - V[0, 1]*V[1, 0])

        AE = 0.25 / area * np.dot(V, V.T)
        BE = area / 12 * (np.ones((3, 3), dtype = np.float64)
                          + np.eye(3, dtype = np.float64))

        for i in range(3):
            for j in range(3):
                A[ele[i], ele[j]] += AE[i, j]
                B[ele[i], ele[j]] += BE[i, j]

    return A, B


# -----------------------
def eigenfunctions(mesh, bnd):

    A, B = fe_matrices(mesh)
    AD = A[bnd == 0, :]
    AD = AD[:, bnd == 0]
    BD = B[bnd == 0, :]
    BD = BD[:, bnd == 0]

    z, V = eigsh(AD, 6, BD, ncv = 100, which = 'SM')

    nn = np.shape(A)[0]
    nnd = np.shape(AD)[0]
    indx = np.zeros(nnd, dtype = np.int32)
    k = 0
    for i in range(nn):
        if bnd[i] == 0:
            indx[k] = i
            k += 1
    
    U = np.zeros((nn, 6), dtype = np.float64)
    for k in range(nnd):
        U[indx[k], :] = V[k, :]

    return z, U
