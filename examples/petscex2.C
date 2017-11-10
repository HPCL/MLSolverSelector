#include "api.h"

/** \file ex1.C
 * \brief Example of using the Petsc Wrappers for the solver selecter and the input file 
 **/

static char help[] = "Solves a tridiagonal linear system with Solver selecter.\n\n";

/** Main */
int main(int argc,char **args)
{
  
  std::string inputfile = "petsc.input";
  std::string database_name = "petsc.database";
  int build_database_from_file = 0;
  int build_database_inline = 0;
  int solver_selecter_on = 1;
  int matrix_dump = 0;

  int ii = 1;
  while ( ii < argc )
  {
      if ( std::string(args[ii]) == "--buildi" )
      {
          ii++;
          build_database_inline = std::stoi(args[ii++]); 
      }
      else if ( std::string(args[ii]) == "--buildf" )
      {
          ii++;
          build_database_from_file = std::stoi(args[ii++]); 
      }
      else if ( std::string(args[ii]) == "--help" )
      {
          std::cout << "This is not much help. \n";
          return 0;    
      }
      else
          std::cout << "Command " << args[ii++] << " not found \n";
  }

  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  Vec            x, b, u;      /* approx solution, RHS, exact solution */
  KSP            ksp;         /* linear solver context */
  PetscReal      norm;  /* norm of solution error */
  PetscScalar    neg_one      = -1.0,one = 1.0,value[3];
 
  std::string filename("petsc.matrices/1138_bus.pbin");
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

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  
  /**************************************************************
   * ******************* Use the solver selecter ***************/
  
  auto interface = std::make_shared<petscUI::PetscUI>(database_name, inputfile, solver_selecter_on,matrix_dump,build_database_inline);

  if ( build_database_from_file ) 
  {
      /* If we are building a database based on matrix files, we simply contruct a solver selecter,
       * Initialize it, and call BuildDataBaseFromFile. */
      petscUI::PetscSolverSelecter ss;
      ss.Initialize(interface);
      ss.BuildDataBaseFromFile(inputfile);
      ss.Finalize();
      ierr = PetscFinalize();
      return 0;
  }  
  
  /* Otherwise, we need to initialize the petsc solver  */
  petscUI::InitializePetscSolverSelecter(ksp,interface);
  
  /* At this point, the ksp is the solver selecter. So, you can continue as 
   * normal. */
  
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);

  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}
