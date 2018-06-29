
#include "Solvers.h"

namespace SolverSelecter
{

Solver::Solver()
{
  solver = "NONE";
  preconditioner = "NONE";
}

Solver::Solver(const std::string &solver,
               const std::string &precond)
{
    SetSolverName(solver,precond);
}

Solver::Solver( const std::string &pstring )
{
    ParseSolverString( pstring );
}


ErrorFlag
Solver::GetSolverInfo(std::string &solver_name,
                      std::string &preconditioner_name,
                      std::set< std::string >  &keys  ) const
{
    solver_name = solver;
    preconditioner_name = preconditioner ;

    keys.clear();
    for ( auto &it : parameters )
        keys.insert(it.first);

    return error_flag;
}

ErrorFlag
Solver::SetSolverName(const std::string &_solver,
                      const std::string &_preconditioner )
{
    if (_solver.size() > 0 )
        solver = _solver;
    else
        solver = "NONE";

    if (_preconditioner.size() > 0 )
        preconditioner = _preconditioner;
    else
        preconditioner = "NONE";
    return error_flag;
}

ErrorFlag
Solver::SetParameter(const std::string &key,
                     const std::string &value )
{
    parameters[key] = value;
    return error_flag;
}

ErrorFlag
Solver::GetParameter(const std::string &key,
                     std::string &value)
{
    auto it = parameters.find(key);
    if (it != parameters.end() )
        value = it->second;
    else
        std::cout << " Parameter key" << key << " not found \n";

    return error_flag;
}

ErrorFlag
Solver::GetSolverString(std::string &pstring ) const 
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

    return error_flag;
}

ErrorFlag
Solver::ParseSolverString( const std::string &pstring )
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

    return error_flag;
}

ErrorFlag
Solver::Print() const 
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
    return error_flag;
}

ErrorFlag
Solver::Clear()
{
    solver = "NONE";
    preconditioner = "NONE";
    parameters.clear();
    return error_flag;
}

}
