
#include <petsc.h>
#include <petscksp.h>
#include <slepceps.h>

#include "mesh.h"
#include "fem.h"

static char help[] = "Salutations.\n\n";

int main(int argc, char **argv)
{
  Mesh m;
  m = ReadMesh(argv[1]);

  int nn = m.nn;
  int *indices = malloc( nn * sizeof(int) );
  for (int n = 0; n < nn; ++n) indices[n] = n;

  SlepcInitialize(&argc, &argv, (char *)0, help);
  int size, me;

  int istart, iend;
  Vec u;
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, nn);
  VecSetFromOptions(u);
  VecGetOwnershipRange(u, &istart, &iend);
  VecAssemblyBegin(u);
  VecAssemblyEnd(u);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &me);

  Mat A = StiffnessMatrix (&m, istart, iend);
  Mat B = MassMatrix      (&m, istart, iend);
  Mat R = NeumannMatrix   (&m, istart, iend);

  int ierr;
  EPS eps; // Eigenvalue problem solver
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, A, NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  int its;
  EPSGetIterationNumber(eps, &its);
  PetscPrintf(PETSC_COMM_WORLD, "Number of iterations of the method: %D\n", its);

  VecDestroy(&u);
  MatDestroy(&A);
  MatDestroy(&B);
  MatDestroy(&R);

  SlepcFinalize();

  free(indices);
  FreeMesh(&m);
}
