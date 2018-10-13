#ifndef SS_SOLVERSELECTER_H
#define SS_SOLVERSELECTER_H

#include "typedefs.h"
#include "MachineLearningInterface.h"
#include "UserInterface.h"
#include "InputFileInterface.h"
#include "DatabaseInterface.h"

namespace SolverSelecter
{

template <typename Matrix, typename Vector>
class SolverSelecter
{

public:
    std::shared_ptr<MachineLearningInterface> machinemodel;
    std::shared_ptr<DatabaseInterface >  database;
    std::shared_ptr<UserInterface<Matrix,Vector>> interface;
    
    std::vector< Solver > solvers;
    std::vector< std::string > matrix_filenames; 
    Solver solver;

    SolverSelecter( std::shared_ptr< UserInterface<Matrix,Vector>> _interface);

    ErrorFlag Initialize(std::string inputfile);

    virtual ~SolverSelecter();

    // Build a database based on an input file. 
    ErrorFlag BuildDataBaseFromFile();

    ErrorFlag ConvertArffFileToDatabase(); 


    // Cross validate a machine learning model for a variety of algorithms 
    ErrorFlag CrossValidate( std::vector< std::string >  algorithm , bool all);

    /** Print the parameter set for the chosen solver */
    ErrorFlag PrintSolver( Solver &solver);

    /** Build the machine learning model and serialize it */
    ErrorFlag BuildModelAndSerialize( std::string serial_name ); 

    /** Solve a system using the solver selecter */
    ErrorFlag Solve(Matrix &A, Vector &x, Vector &b );

    /** Classify the solvers in the database. matrix = -1 signifies all matrices should
     * be re-classified. */
    ErrorFlag ClassifySolvers(int matrix = -1);

private:

    ErrorFlag
    MeasuredSolve(Matrix &A,
                  Vector &x,
                  Vector &b,
                  Solver &solver,
                  std::map<std::string,double> &mresult );

    ErrorFlag
    InitTestSystem(std::string filename,
                   std::map<std::string,double> &fmap,
                   std::unique_ptr< Matrix >    &A,
                   std::unique_ptr< Vector >    &x,
                   std::unique_ptr< Vector >    &b );

    ErrorFlag
    ParameterSpaceFromFile(std::string &database_name,
                           std::vector< std::string > &matrix_filenames,
                           std::vector< Solver > &solvers );

    ErrorFlag
    FreeTestVectors(std::unique_ptr< Vector >  &x,
                    std::unique_ptr< Vector >  &b );

    ErrorFlag
    FreeTestSystem(std::unique_ptr< Matrix >  &A,
                   std::unique_ptr< Vector >  &x,
                   std::unique_ptr< Vector >  &b );
    ErrorFlag
    BuildDataBaseInline(Matrix &A, Vector &x, Vector &b );

};

}

#include "SolverSelecter.tpp"

#endif
