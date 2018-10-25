
#ifndef _SS_FSOLVERS_H
#define _SS_FSOLVERS_H

#include "typedefs.h"

namespace SolverSelecter
{

  /** General struct for defining a solver */
class Solver
{
public:

    std::string preconditioner; /* name of the preconditioner */
    std::string solver; /* Name of the solver */
    std::map< std::string, std::string > parameters; /* parameter key-value pairs */

    /** Constructor */
    Solver( );


    /** Constructor */
    Solver(const std::string &solver,
           const std::string &precond );

    /** Constructor */
    Solver( const std::string &pstring );

    /** Get all the information we have about the solver */
    ErrorFlag
    GetSolverInfo(std::string &solver_name,
                  std::string &preconditioner_name,
                  std::set< std::string >  &keys  ) const;

    /** Change the name of the linear solver */
    ErrorFlag
    SetSolverName(const std::string &_solver,
                  const std::string &_preconditioner );
    /** Change a parameter value */
    ErrorFlag
    SetParameter( const std::string &key,
                  const std::string &value );
    /** Get a parameter value. */
    ErrorFlag
    GetParameter( const std::string &key,
                  std::string &value);

    /** Set a single string that defines the solver -- i.e., serialize the solver to string */
    ErrorFlag
    GetSolverString(std::string &pstring ) const;

    /** Deserialize the solver from string */
    ErrorFlag
    ParseSolverString( const std::string &pstring );

    /** Print the solver information to stdout **/
    ErrorFlag Print() const;

    /** Reset the solve -- in the reset case, the default solver is used. */
    ErrorFlag Clear();


};

}
#endif
