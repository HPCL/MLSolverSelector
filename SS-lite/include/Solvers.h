
#ifndef _SS_FSOLVERS_H
#define _SS_FSOLVERS_H

#include <map>
#include <string>
#include <set>


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
    int
    GetSolverInfo(std::string &solver_name,
                  std::string &preconditioner_name,
                  std::set< std::string >  &keys  ) const;

    /** Change the name of the linear solver */
    int
    SetSolverName(const std::string &_solver,
                  const std::string &_preconditioner );

    /** Change a parameter value */
    int
    SetParameter( const std::string &key,
                  const std::string &value );

    /** Get a parameter value. */
    int
    GetParameter( const std::string &key,
                  std::string &value);

    /** Set a single string that defines the solver -- i.e., serialize the solver to string */
    int
    GetSolverString(std::string &pstring ) const;

    /** Deserialize the solver from string */
    int
    ParseSolverString( const std::string &pstring );

    /** Print the solver information to stdout **/
    int Print() const;

    /** Reset the solve -- in the reset case, the default solver is used. */
    int Clear();

};

#endif
