#include "PetscInterface.h"
#include <memory>

/** \file ex1.C
 * \brief Example of using the Petsc to build a database from file.  
 **/

static char help[] = "Build database from file. .\n\n";

/** Main */
int main(int argc,char **args)
{
    

  std::string matvec = "true";
  int interior = 10;
  int edge = 10;
  int solve = 1;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);
  
  std::shared_ptr<SolverSelecter::PetscUI> ss_interface = std::make_shared< SolverSelecter::PetscUI >();  
  auto ss = std::unique_ptr<SolverSelecter::PetscSS>(new SolverSelecter::PetscSS(ss_interface));

  edge = std::atoi(args[1]);
  interior = std::atoi(args[2]);
  std::string matrix = args[3];  
  solve = std::atoi(args[4]);
  //Set the input parameters for the interface. In this case, all we need to do is pass the interface 
  //the inputfile. 
    
  std::map<std::string, std::string > parameters;
  parameters["inputfile"] = "./moose.input";
  parameters["petsc.internalSample"] = std::to_string(5);
  parameters["petsc.edgeSample"] = std::to_string(5);
  parameters["petsc.matvec"] = matvec;
  parameters["CNAME"] = "Extraction_5_5";

  std::map<std::string, std::string > parameters2;
  parameters2["inputfile"] = "./moose.input";
  parameters2["petsc.internalSample"] = std::to_string(10);
  parameters2["petsc.edgeSample"] = std::to_string(10);
  parameters2["petsc.matvec"] = matvec;
  parameters2["CNAME"] = "Extraction_10_10";
   
  std::map<std::string, std::string > parameters3;
  parameters3["inputfile"] = "./moose.input";
  parameters3["petsc.internalSample"] = std::to_string(10);
  parameters3["petsc.edgeSample"] = std::to_string(5);
  parameters3["petsc.matvec"] = matvec;
  parameters3["CNAME"] = "Extraction_10_5";
   
  std::vector<std::map<std::string,std::string>> pvec;
  pvec.push_back(parameters);
  pvec.push_back(parameters2);
  pvec.push_back(parameters3);

  std::ostringstream oss;
  oss << "./featureAnalysis_" << interior << "_" << edge << "_" << matvec << ".data";

  ss->FeatureExtractionAnalysis(pvec, matrix, oss.str(), solve );
  
  // Clean up some stuff to avoid MPI calls after MPI finalize. 
  ss_interface.reset();
  ss.reset();  
  
  ierr = PetscFinalize();

  return 0;
}  
