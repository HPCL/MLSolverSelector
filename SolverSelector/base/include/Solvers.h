
#ifndef _SS_FSOLVERS_H
#define _SS_FSOLVERS_H

#include "typedefs.h"

namespace SolverSelecter
{

class Solver
{
public:

    std::string preconditioner;
    std::string solver;
    std::map< std::string, std::string > parameters;

    Solver( );

    Solver(const std::string &solver,
           const std::string &precond );

    Solver( const std::string &pstring );

    ErrorFlag
    GetSolverInfo(std::string &solver_name,
                  std::string &preconditioner_name,
                  std::set< std::string >  &keys  ) const;

    ErrorFlag
    SetSolverName(const std::string &_solver,
                  const std::string &_preconditioner );

    ErrorFlag
    SetParameter( const std::string &key,
                  const std::string &value );

    ErrorFlag
    GetParameter( const std::string &key,
                  std::string &value);

    ErrorFlag
    GetSolverString(std::string &pstring ) const;

    ErrorFlag
    ParseSolverString( const std::string &pstring );

    ErrorFlag
    Print() const;

    ErrorFlag
    Clear();


};

}
#endif
