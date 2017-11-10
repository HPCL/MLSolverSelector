
#include "Solvers.h"

_SS_Solver::_SS_Solver( )
{
    SetSolverName( "NONE","NONE" );
}

_SS_Solver::_SS_Solver( const std::string &solver, /**< input, name of the solver (i.e, gmres) */
                        const std::string &precond /**< input, name of the preconditioner (i.e, boomeramg) */
                      )
{
    SetSolverName(solver,precond);
}

_SS_Solver::_SS_Solver( const std::string &pstring /**< input, parameter string obtained with ParseParameterString */)
{
    ParseSolverString( pstring );
}

_SS_ErrorFlag _SS_Solver::GetSolverInfo( std::string &solver_name,         /**< output, solver name */
        std::string &preconditioner_name, /**< output, preconditioner name */
        std::set< std::string >  &keys    /**< output, keys for the parameters */
                                       ) const
{

    solver_name = solver;
    preconditioner_name = preconditioner ;

    keys.clear();
    for ( auto &it : parameters )
        keys.insert(it.first);

    return _SS_error_flag;
}


_SS_ErrorFlag _SS_Solver::SetSolverName( const std::string &_solver, const std::string &_preconditioner )
{
    if (_solver.size() > 0 )
        solver = _solver;
    else
        solver = "NONE";

    if (_preconditioner.size() > 0 )
        preconditioner = _preconditioner;
    else
        preconditioner = "NONE";

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Solver::SetParameter( const std::string &key, const std::string &value )
{
    parameters[key] = value;
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Solver::GetParameter( const std::string &key, std::string &value)
{
    auto it = parameters.find(key);
    if (it != parameters.end() )
        value = it->second;
    else
        std::cout << " Parameter key" << key << " not found \n";

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Solver::GetSolverString( std::string &pstring /**< output, parameter string */) const
{

    pstring.clear();

    std::ostringstream oss;

    /* add the solver and preconditioner */
    std::string sname = solver;

    oss << sname << " " ;
    oss << preconditioner << " " ;

    std::string value;

    /* add all the parameters */
    for ( auto  &it : parameters )
    {
        oss << it.first << " " << it.second << " ";
    }
    pstring = oss.str();

    return _SS_error_flag;
}


_SS_ErrorFlag _SS_Solver::ParseSolverString( const std::string &pstring /**< input, parameter string from GetPar..String */)
{

    Clear(); /* Clear all the parameters */
    std::string value, param, buff, pname, sname;
    std::stringstream ss(pstring);

    // Get solver name;
    ss >> buff;
    sname = buff;
    solver = sname;
    //Get the preconditioner name
    ss >> buff;
    pname = buff;
    preconditioner = pname;

    while ( ss >> buff )
    {
        param = buff ;
        ss >> buff;
        value = buff ;
        parameters.insert(std::make_pair(param,value));
    }

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Solver::Print() const
{

    std::cout << " \t------------------------------------------------------- \n";
    std::cout << " \t-----------RNET SOVER SELECTER SOLVER STATS------------ \n";
    std::cout << " \t------------------------------------------------------- \n";
    std::cout << " \t\tLinear Solver: " << solver << std::endl;
    std::cout << " \t\tPreconditioner: " << preconditioner << std::endl;
    std::cout << " \t\tParameters: \n";

    for ( auto &it  : parameters )
    {
        std::cout << "\t\t\t" << it.first << " : " << it.second << std::endl;
    }
    std::cout << " \t------------------------------------------------------- \n";
    std::cout << " \t------------------------------------------------------- \n";
    return _SS_error_flag;
}
_SS_ErrorFlag _SS_Solver::Clear()
{
    solver = "NONE";
    preconditioner = "NONE";
    parameters.clear();
    return _SS_error_flag;
}

