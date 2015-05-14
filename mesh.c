#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max(a,b) ((a)<(b) ? (b) : (a))

typedef struct {
    int nn, ne, nl;
    double *x, *y;
    int *bnd, *ele, *neigh, *edge, *bnd_edge;
} Mesh;

Mesh ReadMesh(char *meshname) {
    Mesh m;
    int i, n;

    // Read in the nodes
    char *nodename = malloc(strlen(meshname)+6);
    strcpy(nodename,meshname);
    strcat(nodename,".node");

    FILE *nodes;
    nodes = fopen(nodename,"r");
    free(nodename);
    fscanf(nodes, "%d %d %d %d", &m.nn, &i, &i, &i);
    m.x = (double *)malloc( m.nn * sizeof(double) );
    m.y = (double *)malloc( m.nn * sizeof(double) );
    m.bnd = (int *)malloc( m.nn * sizeof(int) );
    for (n = 0; n<m.nn; n++) {
        fscanf(nodes, "%d %lf %lf %d", &i, m.x+n , m.y+n, m.bnd+n);
    }
    fclose(nodes);

    // Read in the elements
    char *elemname = malloc(strlen(meshname)+5);
    strcpy(elemname,meshname);
    strcat(elemname,".ele");

    FILE *elements;
    elements = fopen(elemname,"r");
    free(elemname);

    fscanf(elements, "%d %d %d", &m.ne, &i, &i);
    m.ele = (int *)malloc( 3*m.ne * sizeof(int) );
    for (n=0; n<m.ne; n++) {
        fscanf(elements, "%d %d %d %d", &i, \
            m.ele+3*n, m.ele+3*n+1, m.ele+3*n+2);
        m.ele[ 3*n ] = m.ele[ 3*n ]-1;
        m.ele[3*n+1] = m.ele[3*n+1]-1;
        m.ele[3*n+2] = m.ele[3*n+2]-1;
    }
    fclose(elements);

    // Read in the neighbors
    char *neighname = malloc(strlen(meshname)+7);
    strcpy(neighname,meshname);
    strcat(neighname,".neigh");

    FILE *neighbors;
    neighbors = fopen(neighname,"r");
    free(neighname);
    fscanf(neighbors, "%d %d", &i, &i);
    m.neigh = (int *)malloc( 3*m.ne * sizeof(int) );
    for (n=0; n<m.ne; n++) {
        fscanf(neighbors, "%d %d %d %d", &i, \
            m.neigh+3*n, m.neigh+3*n+1, m.neigh+3*n+2);
        m.neigh[ 3*n ] = max(m.neigh[ 3*n ]-1,-1);
        m.neigh[3*n+1] = max(m.neigh[3*n+1]-1,-1);
        m.neigh[3*n+2] = max(m.neigh[3*n+2]-1,-1);
    }
    fclose(neighbors);

    // Read in the edges
    char *edgename = malloc(strlen(meshname)+6);
    strcpy(edgename,meshname);
    strcat(edgename,".edge");

    FILE *edges;
    edges = fopen(edgename,"r");
    free(edgename);
    fscanf(edges, "%d %d", &m.nl, &i);
    m.edge = (int *)malloc( 2*m.nl * sizeof(int) );
    m.bnd_edge = (int *)malloc( m.nl * sizeof(int) );
    for (n=0; n<m.nl; n++) {
        fscanf(edges, "%d %d %d %d", &i, \
            m.edge+2*n, m.edge+2*n+1, m.bnd_edge+n);
        m.edge[ 2*n ] = m.edge[ 2*n ]-1;
        m.edge[2*n+1] = m.edge[2*n+1]-1;
    }
    fclose(edges);

    // Return the mesh
    return m;
}


void FreeMesh(Mesh *m)
{
  m->nn = 0;
  m->ne = 0;
  m->nl = 0;

  free(m->x);
  free(m->y);
  free(m->bnd);
  free(m->ele);
  free(m->neigh);
  free(m->edge);
  free(m->bnd_edge);
}


/* Determine if a triangle n contains the point w */
int point_in_triangle( Mesh *m, int n, double *w) {
    int i,j,n1,n2;
    double dot_product;
    double dx[2], dw[2];
    for (i=0; i<3; i++) {
        j = (i+1)%3;
        n1 = (*m).ele[3*n+i];
        n2 = (*m).ele[3*n+j];
        dx[0] = -((*m).y[n2] - (*m).y[n1]);
        dx[1] =   (*m).x[n2] - (*m).x[n1];
        dw[0] = w[0] - (*m).x[n1];
        dw[1] = w[1] - (*m).y[n1];
        dot_product = dx[0]*dw[0]+dx[1]*dw[1];
        if (dot_product < 0) {
            return 0;
        }
    }
    return 1;
}


/* Find which triangle contains a point */
int find_point_in_triangle( Mesh *m, double *w) {
    int i, n;

    for (n=0; n<(*m).ne; n++) {
        /* Count how many edges the point lies on  the left side of */
        i = point_in_triangle(m,n,w);
        if (i == 1) {
            return n;
        }
    }
    return -1;
}


/* Given a triangle n, a point W and a displacement dw, determine which
edge, if any, the line segment (W,W+dw) passes through */
int find_edge_crossed( Mesh *m, int n, double *W, double *dw, \
    double *r ) {
    int i,j,k,n1,n2;
    double a,b,c,d,s,t,det;
    double x[2], y[2];

    for (i=0; i<3; i++) {
        j = (i+1)%3;
        k = (i+2)%3;

        n1 = (*m).ele[3*n+j];
        n2 = (*m).ele[3*n+k];
        x[0] = (*m).x[n1];
        y[0] = (*m).y[n1];
        x[1] = (*m).x[n2];
        y[1] = (*m).y[n2];

        /* Set up a linear system for the point of intersection between
        edge i of the current triangle and the line segment between W
        and W+dw */
        a = dw[0];
        b = -(x[1]-x[0]);
        c = dw[1];
        d = -(y[1]-y[0]);

        det = a*d-b*c;
        s = -1.0;
        t = -1.0;
        /* If the determinant of the system is zero, the two lines are
        co-linear */
        if ( det != 0.0) {
            s = ( d*(x[0]-W[0])-b*(y[0]-W[1]) )/det;
            t = ( a*(y[0]-W[1])-c*(x[0]-W[0]) )/det;

            /* If the system has a solution in the unit box, then the
            two line segments intersect */
            if ( 0<s && s<1 && 0<t && t<1) {
                *r = s;
                return i;
            }
        }

    }

    return -1;
}
