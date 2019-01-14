#include "C50Selector.h"
#include <iostream>
#include <fstream>
#include <sstream>

C50Interface::C50Interface(std::string filestem) : m_filestem(filestem) {
  loadModel();
}

C50Interface::~C50Interface() {};

int C50Interface::loadModel() {
 
 if ( !initialized ) { 
    std::cout << " Loading Model with FileStem " << m_filestem << std::endl;
    predictor.reset( new C50Predictor(m_filestem.c_str()));
        
    // Step 1. Get a map of solver hashes to SS solvers 
    std::string solversFile = m_filestem + ".solvers";

    std::ifstream file (solversFile.c_str());
    std::string value;
    int count = 0;
    int solver = 0;

    std::vector<std::string> result;
    if ( !file.is_open() ) {
      std::cout << " <filestem>.solvers does not exist. This file is needed in addition to \n";
      std::cout << " the standard C50 files. This file is used to map the solver hashes in \n";
      std::cout << " the C50 model to the actual petsc parameter sets that were use to do  \n";
      std::cout << " the actual solve. See the examples directory for the correct format \n";
      std::abort();
    }
    while (!file.eof() ) {
      
      result.clear();
      std::getline( file, value); 
      StringSplit(value, ",", result); 
      if ( result.size() == 2 ) {
        solverMap[std::atoi(result[0].c_str())] = result[1];
      }
    } 
    return 1; 
  }
  return 1; 
}

int C50Interface::StringSplit(const std::string &s, const char *delim, std::vector< std::string > &result )
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while ( std::getline( ss,item,delim[0]) )
    {
        if (!item.empty())
            result.push_back( item );
    }
    return 1;
}


Solver C50Interface::Predict( std::map<std::string, double> &features ) {
    
      Solver chosenSolver;
      
      int solverChoice = predictor->findAGoodSolver(features);

      if ( solverChoice < 0 ) { 
        std::cout << "No good solver was found -- Default solver will be used"; 
        chosenSolver.Clear();
      } else { 
         chosenSolver.ParseSolverString(solverMap[solverChoice]);
      }
      return chosenSolver ;
}




