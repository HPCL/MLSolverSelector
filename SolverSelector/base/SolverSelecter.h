#ifndef SS_SOLVERSELECTER_H
#define SS_SOLVERSELECTER_H

#include "typedefs.h"
#include "MachineLearning.h"
#include "UserInterface.h"

/** \file SolverSelecter.h
 * \brief Header file for the solver selecter class
 */

/**
 * The SolverSelecter class is the main user implimented API for the entire project. The
 * user must impliment a derived class of the solver selecter with Matrix and vector set
 * to their problem specific matricies and vector data structures. For example, in petsc,
 * the call  is derived-class : public _SS_SolverSelecter<Mat,Vec>(...). The user must also
 * impliment some pure virtual functions, dicussed below
 * \tparam Matrix The matrix data structure to be used.
 * \tparam Vector The vector data structure to be used.
 **/



template <typename Matrix, typename Vector>
class _SS_SolverSelecter {

public:
    std::shared_ptr<_SS_MachineLearning> machinemodel;                         /**< pointer the the machine learning class */
    std::shared_ptr< _SS_DataBaseBase >  database;                             /**< pointer to the database */
    
    _SS_Solver solver; /** This is the selected solver. Allows user to call back and get the last selected solver */

    std::shared_ptr<_SS_UserInterface<Matrix,Vector>> interface; /**< the user interface implimenting all this stuff */

    /** Constructor */
    _SS_SolverSelecter( std::shared_ptr< _SS_UserInterface<Matrix,Vector>> _interface)
    {
        interface = _interface;
        interface->GetDataBaseImplimentation(database);
        interface->GetMachineLearningModel(machinemodel);
        machinemodel->Initialize(database);
        database->Initialize();

    }

    /** Destructor */
    virtual ~_SS_SolverSelecter()
    {
        if ( interface->build_inline )
            ClassifySolvers();

        database->Finalize();
        database.reset();
        machinemodel.reset();
        interface.reset();
    }

    /** This is a convienece feature that allows the user to add rows to the database using an input file
     * and a list of saved matricies. In this case, the user must supply a function that opens the matrix
     * from file and initializes it. Given this, we can recurse over the parameter space listed in the input
     * file, makeing measurements as we go. ( This is what used to be the "DataBaseBuilder" );
     **/
    _SS_ErrorFlag BuildDataBaseFromFile()
    {

        std::unique_ptr<Matrix> A;
        std::unique_ptr<Vector> x, b;
        std::vector< _SS_Solver > solvers;
        std::vector< std::string > matrix_filenames;
        ParameterSpaceFromFile( interface->input_file, matrix_filenames, solvers );

        for ( auto matrix_file : matrix_filenames )
        {
            std::map<std::string, double> features_map;        
            InitTestSystem( matrix_file, features_map, A, x, b );
            for ( auto next_solver : solvers )
            {
                std::map<std::string, double> measurements_map;
                MeasuredSolve( *A, *x, *b, next_solver, measurements_map );
                database->AddRow( next_solver, matrix_file, features_map, measurements_map );
            }
            FreeTestSystem( A, x, b );
        }

        ClassifySolvers();
        return _SS_error_flag;
    }

    _SS_ErrorFlag CrossValidate( std::vector< std::string >  algorithm , bool all)
    {
        machinemodel->CrossValidateAll( algorithm, all );
    }

    /** Print the parameter set for the chosen solver */
    _SS_ErrorFlag PrintSolver( _SS_Solver &solver) {
        solver.Print();
        return _SS_error_flag;
    }


    /** Select solver and solve the linear system. This is the main public function used
     * to solve the system. Inside here, we do the classification, ml, etc  */
    _SS_ErrorFlag Solve( Matrix &A, Vector &x, Vector &b ) {

        if ( interface->matrix_dump ) 
          interface->DumpMatrix( A );
        
        if ( interface->build_inline)
            BuildDataBaseInline( A , x, b );

        if ( interface->build_inline || ! interface->use_ml )
        {
            interface->GetDefaultSolver(solver) ;
            interface->SolveSystem(A, x, b, solver);
        }
        else
        {
            std::map<std::string, double> features_map;
            machinemodel->Train();
            interface->ExtractFeatures( A, features_map);
            machinemodel->Classify(features_map, solver );

            if (solver.solver == "NONE" )
                interface->GetDefaultSolver(solver);

            interface->SolveSystem( A, x, b, solver );
        }

        return _SS_error_flag;
    }


    /**
     * Classify a solver as good as bad
     **/
    _SS_ErrorFlag ClassifySolvers( ) 
    {
        std::map<std::string, double > class_values; 
        interface->ClassificationMeasurements( class_values );
        database->ClassifySolvers( class_values );
        return _SS_error_flag;
    }

private:

    _SS_ErrorFlag GetSolvers( std::string solver, std::string preconditioner,
                              _SS_parameters_map &sparameters,
                              _SS_parameters_map &pparameters,
                              std::vector< _SS_Solver > &solver_list )
    {
        _SS_parameters_map pmap;
        pmap.insert(sparameters.begin(),sparameters.end());
        pmap.insert(pparameters.begin(),pparameters.end());

        std::vector< _SS_Solver > temp_list;
        if ( pmap.size() > 0 )
        {
            RecurseParameterSpace(solver, preconditioner,pmap, temp_list );
            for (auto it : temp_list )
                solver_list.push_back( it );
            solver_list.pop_back();
        }
        else
            solver_list.push_back( _SS_Solver(solver,preconditioner) );

        return _SS_error_flag;
    }


    _SS_ErrorFlag RecurseParameterSpace(std::string solver, std::string precond,
                                        _SS_parameters_map solver_params,
                                        std::vector< _SS_Solver > &solver_list )
    {


        if (solver_list.size() == 0 )
            solver_list.push_back( _SS_Solver(solver,precond) );

        if ( solver_params.size() == 0 )
        {
            std::string pstring;
            solver_list.back().GetSolverString( pstring ) ;
            solver_list.push_back(_SS_Solver( pstring ));
            return _SS_error_flag;
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

        return _SS_error_flag;
    }

    _SS_ErrorFlag ParameterSpaceFromFile(std::string fname,
                                         std::vector< std::string > &filenames,
                                         std::vector< _SS_Solver >  &solver_list  )
    {

        std::fstream file;
        file.open(fname.c_str(), std::ios::in | std::ios::out | std::ios::app | std::ios::binary );

        std::string line, l;
        std::vector <std::string > split_line, result;

        //adding = -1 for something else, 0 for solver, 1 for precond
        int adding = -1;

        std::string solver_name;
        std::map< std::string, std::set< std::string > > solver_params;
        std::map< std::string, std::map< std::string, std::set< std::string > > > solvers;
        std::map< std::string, std::map< std::string, std::set< std::string > > > precons;
        std::set< std::pair< std::string, std::string > > solver_pairs;

        while (!file.eof())
        {
            std::getline(file,line); /* get the next line in the file */
            result.clear();

            _SS_Utils::StringSplit( line, " ", result );

            if ( result.size() == 0 )
            {
                /* Do nothing */
            }
            else if ( result[0] == "@PARAMETER" )
            {
                if ( adding < 0 )
                    std::cout << " File in wrong format. Adding parameter without solver \n" ;
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
                adding = -1;

                if ( result[0] == "@MATRIX" )
                {
                    filenames.push_back(result[1]);
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
                else
                {
                    std::string fi(result[0].begin(),result[0].begin()+1);
                    if ( fi != "#" )
                        std::cout << " Keyword " << result[0] << " not recongnized \n";
                }
            }
        }
        /* Finally, add all the solvers */
        for ( auto it : solver_pairs )
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

        if ( filenames.size() == 0 || solver_list.size() == 0  )
        {
            std::cout << " Does the file contain a solver and a matrix ???. Also, Error handling? TODO " << std::endl;
        }

        return _SS_error_flag;
    }

    /** Measured solve of the system */
    _SS_ErrorFlag MeasuredSolve( Matrix &A, Vector &x, Vector &b, _SS_Solver &solver, std::map<std::string,double> &mresult ) {
        interface->StartMeasurements(A,x,b);
        interface->SolveSystem(A, x, b, solver);
        interface->StopMeasurements(A,x,b,mresult);
        return _SS_error_flag;
    }

    _SS_ErrorFlag InitTestSystem( std::string                   filename,
                                  std::map<std::string,double> &fmap,
                                  std::unique_ptr< Matrix >    &A,
                                  std::unique_ptr< Vector >    &x,
                                  std::unique_ptr< Vector >    &b )
    {
        interface->InitMatrix( filename, A);
        interface->ExtractFeatures( *A, fmap );
        interface->InitVector( *A, x );
        interface->InitVector( *A, b );
        interface->SetVector( x, "rand" );
        interface->SetVector( b, "ones" );
        return _SS_error_flag;
    }


    _SS_ErrorFlag FreeTestVectors(std::unique_ptr< Vector >  &x,
                                  std::unique_ptr< Vector >  &b )
    {
        interface->FreeVector(x);
        interface->FreeVector(b);
        x.reset(nullptr);
        b.reset(nullptr);
        return _SS_error_flag;
    }


    _SS_ErrorFlag FreeTestSystem( std::unique_ptr< Matrix >  &A,
                                  std::unique_ptr< Vector >  &x,
                                  std::unique_ptr< Vector >  &b )
    {
        interface->FreeMatrix(A);
        FreeTestVectors(x,b);
        A.reset(nullptr);
        return _SS_error_flag;
    }

    _SS_ErrorFlag BuildDataBaseInline( Matrix &A, Vector &x, Vector &b ) {
        
        std::vector< _SS_Solver > solvers;
        std::vector< std::string > matrix_filenames;
        ParameterSpaceFromFile( interface->input_file, matrix_filenames, solvers );
             
        std::string matrix_name;
        std::map<std::string, double> features_map;
        interface->ExtractFeatures( A, features_map);
        interface->GetMatrixName( A, matrix_name );

        std::unique_ptr<Vector> xx, bb;
        interface->CopyVector( x, xx );
        interface->CopyVector( b, bb );
        interface->SetVector( xx, "rand" );
        interface->SetVector( bb, "ones" );

        for ( auto next_solver : solvers  )
        {
            std::map<std::string, double> measurements_map;
            MeasuredSolve( A, *xx, *bb, next_solver, measurements_map );
            database->AddRow( next_solver, matrix_name, features_map, measurements_map );
        }

        FreeTestVectors(xx,bb);
        ClassifySolvers();
        return _SS_error_flag;
    }

};


#endif
