#include "PetscUI.h"

/** \file ex1.C
 * \brief Example of using the Petsc to build a database from file.  
 **/

static char help[] = "Solves a tridiagonal linear system with Solver selecter.\n\n";

/** Main */
int main(int argc,char **args)
{
    
  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  std::string input_file = "../inputs/petsc.input";
  std::string database_name = "../databases/petsc.database";

  std::shared_ptr<PetscUI> ss_interface = std::make_shared< PetscUI >();
  ss_interface->database_name = database_name;
  ss_interface->input_file = input_file;
  PetscSS ss(ss_interface);

  ss.BuildDataBaseFromFile();
  ierr = PetscFinalize();
  return 0;
  }  
