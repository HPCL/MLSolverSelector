#include "PetscInterface.h"

/** \file ex2.C
 * \brief Example of using the Petsc Wrappers for the solver selecter and the input file 
 **/

static char help[] = "Solves a matrix stored to file with Solver selecteri available through \
                      the command line.i I.e, ./ex02.ss -ksp_type KSPSS -help to see options \n \n";


/** Main */
int main(int argc,char **args)
{

  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  Vec            x, b, u;      /* approx solution, RHS, exact solution */
  PetscReal      norm;  /* norm of solution error */
  PetscScalar    neg_one      = -1.0,one = 1.0,value[3];
 
  std::string filename("./matrix/spd-real-int32-float32");
  
  Mat A;
  PetscViewer fd;
  MatCreate(PETSC_COMM_WORLD, &A);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &fd);
  MatSetFromOptions(A);
  MatLoad(A,fd);
  PetscViewerDestroy(&fd);

  PetscInt n,m;
  MatGetSize(A, &m, &n);
  VecCreate(MPI_COMM_WORLD, &x );
  VecSetSizes(x, PETSC_DECIDE, m );
  VecSetFromOptions(x);  
  VecDuplicate(x,&b);
  VecDuplicate(x,&u);

  ierr = VecSet(u,one);CHKERRQ(ierr);
  ierr = MatMult(A,u,b);CHKERRQ(ierr);
   
  /* This is it -- the only change */
  SolverSelecter::PetscCoupler::CreateSolverSelectorKSP();

  KSP ksp; 
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);  
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return 0;
}
