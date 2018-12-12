#include "PetscInterface.h"
#include <memory>

/** \file ex1.C
 * \brief Example of using the interface to serialize a modeliu based on the database.  
 **/

static char help[] = "Builds a machine learning model based on the given database and serializes it..\n\n";

/** Main */
int main(int argc,char **args)
{
    
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  
  
  std::shared_ptr<SolverSelecter::PetscUI> ss_interface = std::make_shared< SolverSelecter::PetscUI >();
  
  auto ss = std::unique_ptr<SolverSelecter::PetscSS>(new SolverSelecter::PetscSS(ss_interface));
  std::map<std::string, std::string> parameters;
  
  // Set the machine learning interface -- in this case C50 (default is Waffles) 
  parameters["mlinterface"] = "C5.0";
  parameters["waffles.algorithm"] = "RandomForest";
  parameters["C50.filestem"] = "./models/Dec11/RS1_petsc_mfree_UFlorida_artemis_p38_dec11_300";

  ss->CrossValidate(parameters, 10);
  
  ss_interface.reset();
  ss.reset();
  
  ierr = PetscFinalize();
  return 0;
}  
