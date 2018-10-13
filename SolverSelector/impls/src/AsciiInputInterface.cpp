
#include "AsciiInputInterface.h"

namespace SolverSelecter
{


AsciiFileParser::AsciiFileParser( )
    : InputFileInterface( )
{
}

ErrorFlag
AsciiFileParser::Parse(std::vector< std::string > &filenames ,
                       std::vector< Solver > &solver_list )
{


    std::fstream file;

    std::cout << " Opening " << inputfile << std::endl;
    file.open(inputfile.c_str(), std::ios::in | std::ios::out | std::ios::app | std::ios::binary );

    std::string line, l;
    std::vector <std::string > split_line, result;

    //adding = -1 for something else, 0 for solver, 1 for precond
    int adding = -1;

    std::string solver_name;
    std::map< std::string, std::set< std::string > > solver_params, interface_p, database_p, machine_p ;
    std::set< std::pair< std::string, std::string > > solver_pairs;
    std::map< std::string, std::map< std::string, std::set< std::string > > > solvers;
    std::map< std::string, std::map< std::string, std::set< std::string > > > precons; 
    
    while (!file.eof())
    {
        std::getline(file,line); /* get the next line in the file */
        result.clear();

        StringSplit( line, " ", result );

        if ( result.size() == 0 )
        {
            /* Do nothing */
        }
        else if ( result[0] == "@PARAMETER" )
        {
            if ( adding < 0 )
                std::cout << " File in wrong format. Adding parameter without first setting a tag \n" ;
            else
            {
                std::string pname = result[1];
                std::set< std::string > pset;
                if ( result.size() <= 2 )
                    pset.insert("none");
                else
                {
                    for ( unsigned int i = 2; i < result.size(); i++ )
                    {
                        std::string fi(result[i].begin(),result[i].begin()+1);
                        if ( fi == "#" ) break;
                        pset.insert(result[i]);
                    }
                }
                solver_params.insert( std::make_pair(pname, pset) );
            }
        }
        else
        {
            /* Finish up the previous solver if required */
            if ( adding == 0 )
                solvers.insert(std::make_pair( solver_name, solver_params ) );
            else if ( adding == 1 )
                precons.insert(std::make_pair( solver_name, solver_params ) );
            else if (adding == 2 ) 
                interface_p = solver_params;
            else if ( adding == 3 )
                machine_p = solver_params;
            else if ( adding == 4 )
                database_p = solver_params;
            adding = -1;

            if ( result[0] == "@MATRIX" )
            {
                filenames.push_back(result[1]);
            }
            else if ( result[0] == "@MATRIXDIR" )
            {
               
#if WITH_BOOST_FILESYSTEM              
                boost::filesystem::directory_iterator end_itr;
                std::string dirc = result[1];
                std::vector< std::string > extensions(result.begin()+2, result.end());
                for ( auto exe : extensions )
                {
                    if ( boost::filesystem::exists(dirc) && boost::filesystem::is_directory(dirc) )
                    {
                        for ( boost::filesystem::directory_iterator itr(dirc); itr != end_itr; ++itr )
                        { 
                            if (boost::filesystem::is_regular_file(itr->path()))
                            {
                                std::string cfile = itr->path().string();
                                if ( boost::filesystem::extension(cfile) == exe )
                                    filenames.push_back(cfile);
                            }
                        }
                    }
                }
#else
                std::cout << " MATRIXDIR option not supported with out boost filesystem support. \
                             Recompile with boost support to use this feature. Ignoring for now \n";
#endif 
            }
            else if ( result[0] == "@PAIR" )
            {
                for ( unsigned int i = 2; i < result.size(); i++ )
                    solver_pairs.insert( std::make_pair( result[1],result[i] ) );
            }
            else if ( result[0] == "@SOLVER" )
            {
                solver_name = result[1];
                solver_params.clear();
                adding = 0;
            }
            else if ( result[0] == "@PRECON" )
            {
                solver_name = result[1];
                solver_params.clear();
                adding = 1;
            }
            else if ( result[0] == "@INTERFACE" )
            {
                solver_params.clear();
                adding = 2;
            }
            else if ( result[0] == "@MACHINELEARNING" )
            {
                solver_params.clear();
                adding = 3;
            }
            else if ( result[0] == "@DATABASE" )
            {
                solver_params.clear();
                adding = 4;
            }
            else
            {
                std::string fi(result[0].begin(),result[0].begin()+1);
                if ( fi != "#" )
                    std::cout << " Keyword " << result[0] << " not recongnized \n";
            }
        }
    }
   
    GetAllSolvers( solver_pairs, solvers, precons, solver_list );
    
    for ( auto it : interface_p ) interface->SetParameter(it.first, *it.second.begin() );
    for ( auto it : machine_p ) machinemodel->SetParameter(it.first, *it.second.begin() );
    for ( auto it : database_p ) database->SetParameter(it.first, *it.second.begin() );


    if ( filenames.size() == 0 || solver_list.size() == 0  )
    {
        std::cout << " Does the file contain a solver and a matrix ???. Also, Error handling? TODO " << std::endl;
    }
    return error_success;

}
ErrorFlag
AsciiFileParser::GetAllSolvers( std::set< std::pair< std::string, std::string > > &solver_pairs,
                                std::map< std::string, std::map< std::string, std::set< std::string > > > &solvers,
                                std::map< std::string, std::map< std::string, std::set< std::string > > > &precons,
                                std::vector< Solver > &solver_list )
{    


    for ( auto it : solver_pairs ) {
      std::cout << " Pairs " << it.first << " " << it.second << std::endl; 
    }
    for ( auto it : solvers ) {
      std::cout << "Solvers " << it.first << std::endl;
    }
    for ( auto it : precons ) {
      std::cout << " Precons " << it.first << std::endl; 
    }

    /* Finally, add all the solvers */
    for ( auto it : solver_pairs)  
    {
        auto search_solver = solvers.find(it.first);
        auto search_precon = precons.find(it.second);





        if ( search_solver != solvers.end() )
        {
            if ( search_precon != precons.end() )
            {
                GetSolvers( it.first, it.second, search_solver->second, search_precon->second, solver_list );
            }
            else
            {
                std::cout << " precon " << it.second << " must be added to input file \n ";
            }
        }
        else
            std::cout << " solver " << it.first << " must be added to input file \n ";
    }

 
    return error_success;
}


ErrorFlag
AsciiFileParser::GetSolvers( std::string solver, std::string preconditioner,
                             parameters_map &sparameters,
                             parameters_map &pparameters,
                             std::vector< Solver > &solver_list )
{
    parameters_map pmap;
    pmap.insert(sparameters.begin(),sparameters.end());
    pmap.insert(pparameters.begin(),pparameters.end());

    std::vector< Solver > temp_list;
    if ( pmap.size() > 0 )
    {
        RecurseParameterSpace(solver, preconditioner,pmap, temp_list );
        for (auto it : temp_list )
            solver_list.push_back( it );
        solver_list.pop_back();
    }
    else
        solver_list.push_back( Solver(solver,preconditioner) );

    return error_success;
}

ErrorFlag AsciiFileParser::RecurseParameterSpace(std::string solver, std::string precond,
        parameters_map solver_params,
        std::vector< Solver > &solver_list )
{


    if (solver_list.size() == 0 )
        solver_list.push_back( Solver(solver,precond) );

    if ( solver_params.size() == 0 )
    {
        std::string pstring;
        solver_list.back().GetSolverString( pstring ) ;
        solver_list.push_back(Solver( pstring ));
        return error_success;
    }

    auto its = std::prev( solver_params.end() );
    if ( its->second.size() > 0 )
    {
        for ( auto it : its->second )
        {
            solver_list.back().SetParameter( its->first, it );
            auto sparams = solver_params;
            sparams.erase(std::prev( sparams.end() ) );
            RecurseParameterSpace(solver,precond,sparams,solver_list);
        }
    }
    else
    {
        auto sparams = solver_params;
        sparams.erase(std::prev(sparams.end()));
        RecurseParameterSpace(solver,precond,sparams,solver_list);
    }

    return error_success;
}

}

