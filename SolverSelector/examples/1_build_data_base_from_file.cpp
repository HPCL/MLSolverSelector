#include "PetscInterface.h"
#include <memory>

/** \file ex1.C
 * \brief Example of using the Petsc to build a database from file.  
 **/

static char help[] = "Build database from file. .\n\n";

/** Main */
int main(int argc,char **args)
{
    
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  
  std::shared_ptr<SolverSelecter::PetscUI> ss_interface = std::make_shared< SolverSelecter::PetscUI >();  
  auto ss = std::unique_ptr<SolverSelecter::PetscSS>(new SolverSelecter::PetscSS(ss_interface));
  
  // Print the help message ( todo -h command line option ) 
  ss_interface->DumpParameterHelpMessage();

  //Set the input parameters for the interface. In this case, all we need to do is pass the interface 
  //the inputfile. 
  std::map<std::string, std::string > parameters;
  parameters["inputfile"] = "./inputfiles/petsc.input";
  
  // Build the database from file. 
  ss->BuildDataBaseFromFile(parameters);    
  
  // Clean up some stuff to avoid MPI calls after MPI finalize. 
  ss_interface.reset();
  ss.reset();  
  
  ierr = PetscFinalize();

  return 0;
}  
