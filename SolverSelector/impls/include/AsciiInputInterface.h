#ifndef ASCII_INPUTINTERFACE_H
#define ASICC_INPUTINTERFACE_H

#include "InputFileInterface.h"

#if WITH_BOOST_FILESYSTEM
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#endif 

namespace SolverSelecter { 

class AsciiFileParser : public InputFileInterface {

  public: 
  
   AsciiFileParser( );


   ErrorFlag 
   Parse(std::vector< std::string > &filenames,
         std::vector< Solver > &solver_list ); 

  
  protected: 
   ErrorFlag GetAllSolvers( std::set< std::pair< std::string, std::string > > &solver_pairs,
                            std::map< std::string, std::map< std::string, std::set< std::string > > > &solvers,
                            std::map< std::string, std::map< std::string, std::set< std::string > > > &precons,
                            std::vector< Solver > &solver_list );

  
  private:
    
  ErrorFlag 
  GetSolvers(std::string solver,
             std::string preconditioner,
             parameters_map &sparameters,
             parameters_map &pparameters,
             std::vector< Solver > &solver_list );
    
   ErrorFlag
   RecurseParameterSpace(std::string solver,
                         std::string precond,
                         parameters_map solver_params,
                         std::vector< Solver > &solver_list );

};

}
#endif
