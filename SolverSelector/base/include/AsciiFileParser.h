#ifndef ASCII_INPUTINTERFACE_H
#define ASICC_INPUTINTERFACE_H

#include "Solvers.h"

#if WITH_BOOST_FILESYSTEM
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#endif 

namespace SolverSelecter { 

/** Class to Parse the input file */
class AsciiFileParser  {

  public: 
   
   /** Constructor */
   AsciiFileParser( );

   /** Parse the input file returning a list of matrix file names and solvers for testing. */
   ErrorFlag 
   Parse(std::string filename,
         std::vector< std::string > &filenames,
         std::vector< Solver > &solver_list);
  
  private: 

   /** populates the solver list from the pairs, solvers and precons found in the file. */
   ErrorFlag GetAllSolvers( std::set< std::pair< std::string, std::string > > &solver_pairs,
                            std::map< std::string, std::map< std::string, std::set< std::string > > > &solvers,
                            std::map< std::string, std::map< std::string, std::set< std::string > > > &precons,
                            std::vector< Solver > &solver_list );
  
  /** Get all the solver combinations */
  ErrorFlag 
  GetSolvers(std::string solver,
             std::string preconditioner,
             parameters_map &sparameters,
             parameters_map &pparameters,
             std::vector< Solver > &solver_list );
  
   /** Recusive helper function for getting every solver combination */  
   ErrorFlag
   RecurseParameterSpace(std::string solver,
                         std::string precond,
                         parameters_map solver_params,
                         std::vector< Solver > &solver_list );

};

}
#endif
