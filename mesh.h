typedef struct {
    int nn, ne, nl;
    double *x, *y;
    int *bnd, *ele, *neigh, *edge, *bnd_edge;
} Mesh;

Mesh ReadMesh(char *meshname);

void FreeMesh(Mesh *m);

int point_in_triangle( Mesh *m, int n, double *w );

int find_point_in_triangle( Mesh *m, double *w );

int find_edge_crossed( Mesh *m, int n, double *W, double *dw, \
    double *r );
