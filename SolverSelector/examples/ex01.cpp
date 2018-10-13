
#include "PetscInterface.h"


/** \file ex1.C
 * \brief Example of using the Petsc Wrappers for the solver selecter and the input file 
 **/

static char help[] = "Solves a tridiagonal linear system with Solver selecter.\n\n";

/** Main */
int main(int argc,char **args)
{
  
  
  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  Vec            x, b, u;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  PetscReal      norm;  /* norm of solution error */
  PetscInt       i,n = 10,col[3],its;
  PetscScalar    neg_one      = -1.0,one = 1.0,value[3];
  
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  
  srand(time(NULL));
  double xxp = -1 ;
  
  value[0] = xxp; value[1] = 2.0; value[2] = xxp;
  for (i=1; i<n-1; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  i    = n - 1; col[0] = n - 2; col[1] = n - 1;
  ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = xxp;
  ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecSet(u,one);CHKERRQ(ierr);
  ierr = MatMult(A,u,b);CHKERRQ(ierr);
  
  SolverSelecter::PetscCoupler::CreateSolverSelectorKSP();
  
  KSP ksp;  
  
  PetscLogDouble space;
  PetscMallocGetCurrentUsage(&space);
  std::cout << "Before KSP" << space << std::endl; 
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
 
 
  PetscMallocGetCurrentUsage(&space);
  std::cout << "KSP Created " << space << std::endl; 
 
 
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  
  PetscMallocGetCurrentUsage(&space);
  std::cout << "Operators set" << space << std::endl; 

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  PetscMallocGetCurrentUsage(&space);
  std::cout << "After Solve" << space << std::endl; 

  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  PetscMallocGetCurrentUsage(&space);
  std::cout << "After KSP Destroy" << space << std::endl; 
  
  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  
  ierr = PetscFinalize();
  return 0;
}
