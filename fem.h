
#ifndef FEM_H
#define FEM_H

Mat StiffnessMatrix( Mesh *m, int istart, int iend );

Mat MassMatrix( Mesh *m, int istart, int iend );

Mat NeumannMatrix( Mesh *m, int istart, int iend );

#endif
