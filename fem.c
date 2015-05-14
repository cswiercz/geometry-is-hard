#include <stdio.h>
#include <math.h>
#include <petsc.h>
#include <mpi.h>

#include "mesh.h"


////////////////////////////////////////////////////////////////////////////
// Construct the stiffness matrix                                         //
////////////////////////////////////////////////////////////////////////////
Mat StiffnessMatrix( Mesh *m, int istart, int iend ) {
    // Get the number of nodes and elements of the mesh
    int nn = (*m).nn;     int ne = (*m).ne;     int nl = (*m).nl;

    // Declare some local variables
    int i, j, k, n, ele[3], edge[2];
    double area, D[6], AE[9];


    // PETSc stuff to set up the matrix.
    // This means finding out how many non-zero entries are in each row
    // which correspond to local and non-local unknowns, then pre-allocating
    // the matrix accordingly.
    Mat A;

    int d_nnz[iend-istart], o_nnz[iend-istart], edge_is_local[2];
    for (n=0; n<iend-istart; n++) {
        d_nnz[n] = 1;
        o_nnz[n] = 0;
    }
    for (n=0; n<nl; n++) {
        edge[0] = (*m).edge[ 2*n ];
        edge[1] = (*m).edge[2*n+1];

        edge_is_local[0] = istart <= edge[0] && edge[0] < iend;
        edge_is_local[1] = istart <= edge[1] && edge[1] < iend;

        if ( edge_is_local[0] && edge_is_local[1] ) {
            d_nnz[ edge[0]-istart ] = d_nnz[ edge[0]-istart ]+1;
            d_nnz[ edge[1]-istart ] = d_nnz[ edge[1]-istart ]+1;
        } else if ( edge_is_local[0] ) {
            o_nnz[ edge[0]-istart ] = o_nnz[ edge[0]-istart ]+1;
        } else if ( edge_is_local[1] ) {
            o_nnz[ edge[1]-istart ] = o_nnz[ edge[1]-istart ]+1;
        }
    }

    MatCreateAIJ(PETSC_COMM_WORLD,iend-istart,iend-istart,nn,nn, \
        0,d_nnz,0,o_nnz,&A);

    // Tell PETSc that A is SPD
    MatSetOption(A,MAT_SPD,PETSC_TRUE);


    // Fill in the entries of the stiffness matrix.
    // This means constructing the element stiffness matrix, then adding the
    // element matrix into the global matrix.
    for (n=0; n<ne; n++) {
        // Get the nodes for the current element
        ele[0] = (*m).ele[ 3*n ];
        ele[1] = (*m).ele[3*n+1];
        ele[2] = (*m).ele[3*n+2];

        // Construct the element stiffness matrix
        /*  D = [ x_3-x_2 , y_3-y_2 ;
                  x_1-x_3 , y_1-y_3 ;
                  x_2-x_1 , y_2-y_1 ];
            element stiffness matrix = area/4*D*transpose(D)
            See Grossmann's book for derivation of this formula */
        if ( istart <= ele[0] && ele[0] < iend ) {        
            for (i=0; i<3; i++) {
                j = (i+1)%3;
                k = (i+2)%3;
                D[ 2*i ] = (*m).x[ele[k]]-(*m).x[ele[j]];
                D[2*i+1] = (*m).y[ele[k]]-(*m).y[ele[j]];
            }
            area = 0.5*fabs(D[0]*D[3]-D[1]*D[2]);

            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    AE[3*i+j] = 0.0;
                    for (k=0; k<2; k++) {
                        AE[3*i+j] = AE[3*i+j]+0.25*D[2*i+k]*D[2*j+k]/area;
                    }
                }
            }

            // Add the element stiffness matrix to the global matrix
            MatSetValues(A,3,ele,3,ele,AE,ADD_VALUES);
        }
    }

    // Some PETSc stuff to finalize assembly of the matrix
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    // Return the matrix
    return A;
}



////////////////////////////////////////////////////////////////////////////
// Construct the mass matrix                                              //
////////////////////////////////////////////////////////////////////////////
Mat MassMatrix( Mesh *m, int istart, int iend ) {
    // Get the number of nodes and elements of the mesh
    int nn = (*m).nn;     int ne = (*m).ne;     int nl = (*m).nl;

    // Declare some local variables
    int i, j, k, n, ele[3], edge[2];
    double area, D[6], BE[9];

    // PETSc stuff to set up the matrix.
    // This means finding out how many non-zero entries are in each row
    // which correspond to local and non-local unknowns, then pre-allocating
    // the matrix accordingly.
    Mat B;

    int d_nnz[iend-istart], o_nnz[iend-istart], edge_is_local[2];
    for (n=0; n<iend-istart; n++) {
        d_nnz[n] = 1;
        o_nnz[n] = 0;
    }
    for (n=0; n<nl; n++) {
        edge[0] = (*m).edge[ 2*n ];
        edge[1] = (*m).edge[2*n+1];

        edge_is_local[0] = istart <= edge[0] && edge[0] < iend;
        edge_is_local[1] = istart <= edge[1] && edge[1] < iend;

        if ( edge_is_local[0] && edge_is_local[1] ) {
            d_nnz[ edge[0]-istart ] = d_nnz[ edge[0]-istart ]+1;
            d_nnz[ edge[1]-istart ] = d_nnz[ edge[1]-istart ]+1;
        } else if ( edge_is_local[0] ) {
            o_nnz[ edge[0]-istart ] = o_nnz[ edge[0]-istart ]+1;
        } else if ( edge_is_local[1] ) {
            o_nnz[ edge[1]-istart ] = o_nnz[ edge[1]-istart ]+1;
        }
    }

    MatCreateAIJ(PETSC_COMM_WORLD,iend-istart,iend-istart,nn,nn, \
        0,d_nnz,0,o_nnz,&B);

    // Tell PETSc that the matrix is SPD
    MatSetOption(B,MAT_SPD,PETSC_TRUE);


    // Fill in the entries of the mass matrix.
    // This means constructing the element mass matrix, then adding the
    // element matrix into the global matrix.
    for (n=0; n<ne; n++) {
        // Get the nodes for the current element
        ele[0] = (*m).ele[ 3*n ];
        ele[1] = (*m).ele[3*n+1];
        ele[2] = (*m).ele[3*n+2];

        if (istart <= ele[0] && ele[0] < iend) {
            // Compute the area of the current element
            for (i=0; i<3; i++) {
                j = (i+1)%3;
                k = (i+2)%3;
                D[ 2*i ] = (*m).x[ele[k]]-(*m).x[ele[j]];
                D[2*i+1] = (*m).y[ele[k]]-(*m).y[ele[j]];
            }
            area = 0.5*fabs(D[0]*D[3]-D[1]*D[2]);

            // Construct the element mass matrix
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    BE[3*i+j] = area/12.0;
                }
                BE[4*i] = BE[4*i]+area/12.0;
            }

            // Add the element mass matrix to the global matrix
            MatSetValues(B,3,ele,3,ele,BE,ADD_VALUES);
        }
    }

    // Some PETSc stuff to finalize assembly of the matrix
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

    // Return the matrix
    return B;
}



////////////////////////////////////////////////////////////////////////////
// Construct the neumann Matrix                                           //
////////////////////////////////////////////////////////////////////////////
Mat NeumannMatrix( Mesh *m, int istart, int iend ) {
    // Get the number of nodes and elements of the mesh
    int nn = (*m).nn;    int nl = (*m).nl;

    // Declare some local variables
    int n, edge[2];
    double dx, RE[4];

    // PETSc stuff to set up the matrix.
    // This means finding out how many non-zero entries are in each row
    // which correspond to local and non-local unknowns, then pre-allocating
    // the matrix accordingly.
    Mat R;

    int d_nnz[iend-istart], o_nnz[iend-istart], edge_is_local[2];
    for (n=0; n<iend-istart; n++) {
        d_nnz[n] = (*m).bnd[n+istart];
        o_nnz[n] = 0;
    }
    for (n=0; n<nl; n++) {
        if ( (*m).bnd_edge[n] ) {
            edge[0] = (*m).edge[ 2*n ];
            edge[1] = (*m).edge[2*n+1];

            edge_is_local[0] = istart <= edge[0] && edge[0] < iend;
            edge_is_local[1] = istart <= edge[1] && edge[1] < iend;

            if ( edge_is_local[0] && edge_is_local[1] ) {
                d_nnz[ edge[0]-istart ] = d_nnz[ edge[0]-istart ]+1;
                d_nnz[ edge[1]-istart ] = d_nnz[ edge[1]-istart ]+1;
            } else if ( edge_is_local[0] ) {
                o_nnz[ edge[0]-istart ] = o_nnz[ edge[0]-istart ]+1;
            } else if ( edge_is_local[1] ) {
                o_nnz[ edge[1]-istart ] = o_nnz[ edge[1]-istart ]+1;
            }
        }
    }

    MatCreateAIJ(PETSC_COMM_WORLD,iend-istart,iend-istart,nn,nn, \
        0,d_nnz,0,o_nnz,&R);

    // Tell PETSc that the matrix is SPD
    MatSetOption(R,MAT_SPD,PETSC_TRUE);


    // Fill in the entries of the Neumann matrix
    for (n=0; n<nl; n++) {
        if ( (*m).bnd_edge[n]) {
            // Get the nodes for the current edge
            edge[0] = (*m).edge[ 2*n ];
            edge[1] = (*m).edge[2*n+1];

            if (istart <= edge[0] && edge[0] < iend) {
                // Compute the length of the current edge
                dx = sqrt( pow( (*m).x[edge[1]]-(*m).x[edge[0]],2) \
                          +pow( (*m).y[edge[1]]-(*m).y[edge[0]],2) );

                // Compute this edge's Neumann matrix
                RE[0] = dx/3;   RE[1] = dx/6;
                RE[2] = dx/6;   RE[3] = dx/3;

                // Add the edge matrix to the global matrix
                MatSetValues(R,2,edge,2,edge,RE,ADD_VALUES);
            }
        }
    }

    // Some PETSc stuff to finalize assembly of the matrix
    MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY);

    // Return the matrix
    return R;
}
