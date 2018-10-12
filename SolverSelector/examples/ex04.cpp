#include "PetscInterface.h"
#include <memory>

/** \file ex1.C
 * \brief Example of using ML interface to build then serialize a model.  
 **/

static char help[] = "Solves a tridiagonal linear system with Solver selecter.\n\n";

/** Main */
int main(int argc,char **args)
{
    
  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  MPI_Comm_set_errhandler(PETSC_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  std::string input_file = "petsc.input";
  
  std::shared_ptr<SolverSelecter::PetscUI> ss_interface = std::make_shared< SolverSelecter::PetscUI >();
  
  auto ss = std::unique_ptr<SolverSelecter::PetscSS>(new SolverSelecter::PetscSS(ss_interface));
  ss->Initialize(input_file);
  ss->BuildModelAndSerialize("model.serialized");    
 
  //std::vector<std::string> alg;
  //alg.push_back("DecisionTree");
  //ss->CrossValidate(alg, false); 
  /* Destuctors could handle these; but, there is a MPI call is a destructor somewhere which causes
   * an error at exit. Destructors really shouldn't be making MPI calls, so, we need to figure out
   * whats goin on there. */
  ss_interface.reset();
  ss.reset();
  
  ierr = PetscFinalize();
  return 0;
}  
