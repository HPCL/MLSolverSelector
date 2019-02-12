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
  std::vector <std::string> CNames;
  CNames.push_back("T20Time"); 
  // Set the machine learning interface -- in this case C50 (default is Waffles) 
  //parameters["mlinterface"] = "C50";
  ss->SerializeMachineLearningModel(CNames,parameters, "serializedModel_T20trial");  
  ss_interface.reset();
  ss.reset();
  
  ierr = PetscFinalize();
  return 0;
}  
