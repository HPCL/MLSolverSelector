#ifndef SS_SOLVERSELECTER_H
#define SS_SOLVERSELECTER_H

#include "typedefs.h"
#include "MachineLearningInterface.h"
#include "UserInterface.h"
#include "AsciiFileParser.h"
#include "DatabaseInterface.h"

namespace SolverSelecter
{

template <typename Matrix, typename Vector>

  
/** Main API class for calling the various solverselector functions. */  
class SolverSelecter
{

public:
    std::shared_ptr<MachineLearningInterface> machinemodel;
    std::shared_ptr<DatabaseInterface >  database;
    std::shared_ptr<UserInterface<Matrix,Vector>> interface;
    
    std::vector< Solver > solvers;
    std::vector< std::string > matrix_filenames; 
    Solver solver;
    bool inputParsed = false;
    bool initialized = false;

    /** Constructor */
    SolverSelecter( std::shared_ptr< UserInterface<Matrix,Vector>> _interface);
     
    /** Desctructor */
    virtual ~SolverSelecter();

    /** Build a database based on an input file. */ 
    ErrorFlag BuildDataBaseFromFile(std::map<std::string,std::string> &parameters);

    /** Cross validate a machine learning model for a variety of algorithms. This one
     * does cross validation on   */
    ErrorFlag CrossValidate(std::map<std::string, std::string> &parameters, int folds);

    /** Print the parameter set for the chosen solver */
    ErrorFlag PrintSolver( Solver &solver);

    /** Build the machine learning model and serialize it */
    ErrorFlag SerializeMachineLearningModel(std::map<std::string,std::string> &parameters, std::string output ); 

    /** Solve a system using the solver selecter */
    ErrorFlag Solve(Matrix &A, Vector &x, Vector &b );

    /** Time Feature Extraction for the Matrices in the input file, then use the default solver to compare
     * the time to solve the matrix to the time to solve with the default solver. ( Testing --  */
    ErrorFlag FeatureExtractionAnalysis(std::vector<std::map<std::string,std::string>> &parameters, std::string outputfile, bool solve);
    
    ErrorFlag FeatureExtractionAnalysis( std::vector<std::map<std::string,std::string>> &parameters,
                                                          std::string matrixfile,
                                                          std::string outputfile,
                                                          bool solve);


    /** Initialize the solver selector. This is called by all functions prior to use, but 
     * only done once per instance of the class. */
    ErrorFlag Initialize(std::map<std::string, std::string> &parameters);

private:

    ErrorFlag ParseInputFile();

    /** Perform a solve with Performance measurements turned on. This is used to record things
     * like CPU time for putting into the database. */
    ErrorFlag
    MeasuredSolve(Matrix &A,
                  Vector &x,
                  Vector &b,
                  Solver &solver,
                  std::map<std::string,double> &mresult );

    /** Init the testing system for building database entries */
    ErrorFlag
    InitTestSystem(std::string filename,
                   std::map<std::string,double> &fmap,
                   std::unique_ptr< Matrix >    &A,
                   std::unique_ptr< Vector >    &x,
                   std::unique_ptr< Vector >    &b );

    /** Build the parameter space for database building from the input file */
    ErrorFlag
    ParameterSpaceFromFile(std::string &database_name,
                           std::vector< std::string > &matrix_filenames,
                           std::vector< Solver > &solvers );

    /** Free the test vectors that were build for testing */
    ErrorFlag
    FreeTestVectors(std::unique_ptr< Vector >  &x,
                    std::unique_ptr< Vector >  &b );

    /** Free the testing system that was built for database building */
    ErrorFlag
    FreeTestSystem(std::unique_ptr< Matrix >  &A,
                   std::unique_ptr< Vector >  &x,
                   std::unique_ptr< Vector >  &b );
    
    /** Build database in line. This is called in solve when the option is set to build 
     * database entries inside existing simulations */
    ErrorFlag
    BuildDataBaseInline(Matrix &A, Vector &x, Vector &b );

    /** Classify all the solvers in the database (or just one matrix if index given */
    ErrorFlag ClassifySolvers(int matrix = -1);
};

}

#include "SolverSelecter.tpp"

#endif
