
#ifndef _SS_SOLVERS_H
#define _SS_SOLVERS_H

#include "ParameterBase.h"

/**
 * A solver is just a key value pair. The second parameter is not used.
 **/
class _SS_Solver
{

protected:
    std::string preconditioner;                      /**< name of the preconditioner. The class "name" is the solver name */
    std::string solver;                              /**< string representing the solver name */
    std::map< std::string, std::string > parameters; /**< all the parameters */

public:

    _SS_Solver( ) ;

    /** Constructor with names **/
    _SS_Solver( const std::string &solver, /**< input, name of the solver (i.e, gmres) */
                const std::string &precond /**< input, name of the preconditioner (i.e, boomeramg) */ );

    /**
     * Construct based on string obtained previously with ParseParameterString
     **/
    _SS_Solver( const std::string &pstring /**< input, parameter string obtained with ParseParameterString */);

    /**
     * Get the solver, preconditioner, and parameter keys
     **/
    _SS_ErrorFlag GetSolverInfo( std::string &solver_name,         /**< output, solver name */
                                 std::string &preconditioner_name, /**< output, preconditioner name */
                                 std::set< std::string >  &keys    /**< output, keys for the parameters */ ) const;

    _SS_ErrorFlag SetSolverName( const std::string &_solver, const std::string &_preconditioner );

    /** Set the parameter */
    _SS_ErrorFlag SetParameter( const std::string &key, const std::string &value );

    /** Get the parameter based on key */
    _SS_ErrorFlag GetParameter( const std::string &key, std::string &value);

    /**
     * Get parameter strings for this solver. This builds a string from which the current
     * solver can be rebuilt. ( That string is stored in the database )
     **/
    _SS_ErrorFlag GetSolverString( std::string &pstring /**< output, parameter string */) const;


    /**
     * Parse a solver string obtained at some other point with GetSolverString
     */
    _SS_ErrorFlag ParseSolverString( const std::string &pstring /**< input, parameter string from GetPar..String */);
    /**
     * Print out the solver information
     **/
    _SS_ErrorFlag Print() const;

    /**
     * Reset this solver
     */
    _SS_ErrorFlag Clear();


};

#endif
