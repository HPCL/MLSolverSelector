#ifndef _SS_INPUTFILESPEC_HEADER
#define _SS_INPUTFILESPEC_HEADER

#include "Solvers.h"

namespace SolverSelecter { 

  class InputFileInterface {

    public: 

    std::string inputfile;
    std::shared_ptr< SSBase > interface, machinemodel, database ;
    
    InputFileInterface() {};
    
    virtual ~InputFileInterface() {};

    virtual ErrorFlag Initialize(std::string input, 
                                 std::shared_ptr< SSBase > inter,  
                                 std::shared_ptr< SSBase > machine,
                                 std::shared_ptr< SSBase > data ) 
    {
        interface = inter;
        machinemodel = machine;
        database =data;
        inputfile = input;
        return error_flag;
    }  

    virtual ErrorFlag 
    Parse(std::vector< std::string > &matrices,
          std::vector< Solver > &solver_list ) = 0;

  };
}

#endif



