
//Template class impl to be included at the bottom of the header file. 

namespace SolverSelecter
{

template<typename Matrix, typename Vector>
SolverSelecter<Matrix,Vector>::SolverSelecter( std::shared_ptr< UserInterface<Matrix,Vector>> _interface)
{
    interface = _interface;
}

template<typename Matrix, typename Vector>
ErrorFlag 
SolverSelecter<Matrix,Vector>::Initialize(std::map<std::string, std::string > &parameters) {

  if ( !initialized ) {  
      
    interface->SetParameters(parameters);  
    interface->GetDataBaseImplimentation(database);
    interface->GetMachineLearningModel(machinemodel);
    database->Initialize();
    machinemodel->Initialize(database);
  }
  return  error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::ParseInputFile() {

  if ( !inputParsed) {
    AsciiFileParser parser;
    parser.Parse( interface->GetInputFileName(), matrix_filenames, solvers); 
    inputParsed = true; 
  }
  return 0;
}

template<typename Matrix, typename Vector>
ErrorFlag 
SolverSelecter<Matrix,Vector>::SerializeMachineLearningModel(std::map<std::string,std::string> &parameters, std::string output) {
   
    Initialize(parameters);
    machinemodel->Serialize(output);
    return 0;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::BuildDataBaseFromFile(std::map<std::string,std::string> &parameters )
{
    Initialize(parameters);
    ParseInputFile(); 
    
    std::unique_ptr<Matrix> A;
    std::unique_ptr<Vector> x, b;
    
    for ( auto matrix_file : matrix_filenames )
    {
        int row_hash;
        std::cout << " Starting Matrix " << matrix_file << std::endl;
        std::map<std::string, double> features_map;
        InitTestSystem( matrix_file, features_map, A, x, b );
        for ( auto next_solver : solvers )
        {
            std::map<std::string, double> measurements_map;
            MeasuredSolve( *A, *x, *b, next_solver, measurements_map );
            row_hash = database->AddRow( next_solver, features_map, measurements_map );
        }
        FreeTestSystem( A, x, b );
        ClassifySolvers(row_hash);
    }

    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::CrossValidate(std::map<std::string,std::string> &parameters, int folds)
{
    Initialize(parameters);
    machinemodel->CrossValidate( folds );
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::PrintSolver( Solver &solver)
{

    solver.Print();
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::Solve( Matrix &A,
                       Vector &x,
                       Vector &b )
{

    bool success;
    if ( interface->GetDumpMatrix() )
        interface->DumpMatrix( A );

    if ( interface->GetBuildInline() )
        BuildDataBaseInline( A , x, b );

    bool change_solver;
    interface->ChangeSolver(A, x, b , change_solver );
    if ( interface->GetBuildInline() || ! interface->GetUseML() )
    {
        interface->GetDefaultSolver(solver) ;
        interface->SolveSystem(A, x, b, solver, success);
    }
    else if ( ! change_solver )
    {
        std::cout << "Using the same solver as last time \n" ;
        interface->SolveSystem(A,x,b,solver,success); // Use the last solver
    }
    else
    {

        std::map<std::string, double> features_map;
        std::cout <<  "Extracting Features \n "  << machinemodel << std::endl;;
        interface->ExtractFeatures( A, features_map);
        std::cout <<  "Extracted Features \n "  << machinemodel << std::endl;;
        machinemodel->Classify(features_map, solver );
        std::cout << " Classified Solver \n "; 
        if (solver.solver != "NONE" )
        {
            solver.Print();
            interface->SolveSystem( A, x, b, solver, success );
            int count = 0;
            while (!success && count < 5 )
            {
                count++;
                printf("Yikes, the solver we picked diverged -- adding to banned list and restarting-rerunning with default solver\n");
                machinemodel->AddToBanedListAndTryAgain(features_map,solver);
                if (solver.solver == "NONE" ) {
                    printf("Yikes, this this model sucks rerunning with default solver\n");
                    interface->GetDefaultSolver(solver);
                    interface->SolveSystem(A,x,b,solver,success);
                    if ( !success ) {
                        printf("Even the default failed - Bailing \n" );
                        return -1;
                    }
                } else 
                    interface->SolveSystem(A,x,b,solver,success);

            }
            if (!success) {
                interface->SolveSystem(A,x,b,solver,success);
            }
        }
        else
        {
            std::cout << " Solver selection failed -- Using default \n" ;
            interface->GetDefaultSolver(solver);
            interface->SolveSystem(A,x,b,solver,success);
        }
    }
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::ClassifySolvers(int matrix)
{
    std::map<std::string, double > class_values;
    interface->ClassificationMeasurements( class_values );
    database->ClassifySolvers( class_values, matrix );
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::MeasuredSolve(Matrix &A,
                              Vector &x,
                              Vector &b,
                              Solver &solver,
                              std::map<std::string,double> &mresult )
{

    bool success;
    solver.Print();
    interface->StartMeasurements(A,x,b);
    interface->SolveSystem(A, x, b, solver, success);
    interface->StopMeasurements(A,x,b,mresult);
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::InitTestSystem( std::string                  filename,
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
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::FreeTestVectors(std::unique_ptr< Vector >  &x,
                                std::unique_ptr< Vector >  &b )
{

    interface->FreeVector(x);
    interface->FreeVector(b);
    x.reset(nullptr);
    b.reset(nullptr);
    return error_flag;
}


template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::FreeTestSystem( std::unique_ptr< Matrix >  &A,
                                std::unique_ptr< Vector >  &x,
                                std::unique_ptr< Vector >  &b )
{

    interface->FreeMatrix(A);
    FreeTestVectors(x,b);
    A.reset(nullptr);
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
SolverSelecter<Matrix,Vector>::BuildDataBaseInline( Matrix &A,
                                     Vector &x,
                                     Vector &b )
{
    ParseInputFile();
    
    std::string matrix_name;
    std::map<std::string, double> features_map;
    interface->ExtractFeatures( A, features_map);
    interface->GetMatrixName( A, matrix_name );

    std::unique_ptr<Vector> xx, bb;
    interface->CopyVector( x, xx );
    interface->CopyVector( b, bb );

    int row_hash = -1;
    for ( auto next_solver : solvers  )
    {
        interface->SetVector( xx, "rand" );
        std::map<std::string, double> measurements_map;
        MeasuredSolve( A, *xx, *bb, next_solver, measurements_map );
        row_hash = database->AddRow( next_solver, features_map, measurements_map );
    }
    FreeTestVectors(xx,bb);
    ClassifySolvers(row_hash);
    return error_flag;
}


template<typename Matrix, typename Vector>
SolverSelecter<Matrix,Vector>::~SolverSelecter()
{
    database->Finalize();
    database.reset();
    machinemodel.reset();
    interface.reset();
}

}
