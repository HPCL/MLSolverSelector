
#include "petscksp.h"                                                                                         
#include "petsc/private/kspimpl.h" 
  
static char help[] = "Scan petsc for solver options.\n\n";

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  int size;
  PetscInitialize(&argc,&args,(char*)0,help);
  PETSC_COMM_WORLD = MPI_COMM_WORLD; 
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  KSP ksp;
  KSPCreate(MPI_COMM_WORLD,&ksp);
  KSPSetFromOptions(ksp);
 
  ierr = PetscFinalize();
  return 0;
}
