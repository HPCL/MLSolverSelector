#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Solvers.h"
#include "MachineLearningInterface.h"
#include "DatabaseInterface.h"

namespace SolverSelecter { 
/** Base class for a user interface between the solver package and solver selector */
template<typename Matrix,typename Vector>
class UserInterface : public SSBase
{
public:
   
    /* Constructor -- add some important parameters */
    UserInterface() : SSBase("UserInterface") 
    {
      AddParameter("matrix_dump","false","Dump matrix");
      AddParameter("use_ml","true","Use MachineLearning");
      AddParameter("build_inline","false","Build Database inline");
      AddParameter("inputfile", "solver.input", "Input file for defining solvers");
    } ;
     
    /** Some simple get function calls for convienence */
    bool GetDumpMatrix() { return GetParameterAsBool("matrix_dump") ; } 
    bool GetUseML() { return GetParameterAsBool("use_ml"); }
    bool GetBuildInline() { return GetParameterAsBool("build_inline"); }
    std::string GetInputFileName() { return GetParameter("inputfile"); } 

    /** Set the parameters based on a list */
    ErrorFlag SetParameters(std::map<std::string, std::string> &parameters) { 
       for ( auto &it : parameters ) {
          SetParameter(it.first, it.second); 
       }
       return 0;
    }

    /** This function should solve the system Ax = b using the given solver and set 
     * the success flag accordingly. */
    virtual ErrorFlag
    SolveSystem(Matrix &A,
                Vector &x, 
                Vector &b, 
                Solver &solver, 
                bool &success ) = 0 ;

    /** User should implement this function to look at the system and determine if we should
     * use machine learning to change the solver or if we should use the same solver as last time. 
     * This can be used to do things like only changing solver at the start of a newton solve or
     * on every tenth newton iteration */
    virtual ErrorFlag
    ChangeSolver(Matrix &A,
                 Vector &x,
                 Vector &b ,
                 bool &change_it ) = 0;

    /** This is the function that must extract the features from the matrix. TODO -- tell 
     * the user what features we actually need */
    virtual ErrorFlag
    ExtractFeatures(Matrix &A,
                    std::map<std::string, double> &fmap ) = 0;
   
    /** This function should start any measurements that the user wants to make during the next
     * call to solve ( i.e., start a CPUTimer or power monitor */
    virtual ErrorFlag
    StartMeasurements(Matrix &A,
                      Vector &x,
                      Vector &b ) = 0;
    /** Stop all measurements and return the values in the map. */
    virtual ErrorFlag
    StopMeasurements(Matrix &A,
                     Vector &x,
                     Vector &b,
                     std::map<std::string, double> &mmap ) = 0;
    
    /** Return the database interface you want to use ( currently only sqlite3 is supported ) */
    virtual ErrorFlag GetDataBaseImplimentation(std::shared_ptr< DatabaseInterface > &database) = 0;
    
    /** Return the machine learning interface to use ( currently only Waffles and C5.0 supported ) */
    virtual ErrorFlag GetMachineLearningModel(std::shared_ptr< MachineLearningInterface > &machine) = 0;
    
    /** Return a list of bad values to use to classify the solvers. */
    virtual ErrorFlag ClassificationMeasurements( std::map<std::string, double> &class_values ) ; 
    
    /** Implement a function to dump matrix to file */
    virtual ErrorFlag DumpMatrix(Matrix &A ) ;

    /** IMplement function to initialize a matrix from a file on disk (optional)*/ 
    virtual ErrorFlag InitMatrix(std::string filename, std::unique_ptr<Matrix> &A );

    // Implement a function to initialize a vector that can be multiplied against the matrix A;
    virtual ErrorFlag InitVector(const Matrix &AA, std::unique_ptr<Vector> &x );
    
    /** Copy vector AA into x; */
    virtual ErrorFlag CopyVector(const Vector &AA, std::unique_ptr<Vector> &x );
    
    /** Set the values of cloned_x based on type_ TODO */
    virtual ErrorFlag SetVector(std::unique_ptr<Vector> &cloned_x, std::string type_ ) ;

    /** Free the matrix */
    virtual ErrorFlag FreeMatrix( std::unique_ptr<Matrix> &A ) ;

    /** Free the igiven vector */
    virtual ErrorFlag FreeVector( std::unique_ptr<Vector> &x ) ;
    
    /** Return a name for the matrix. This is used to set the filename when dumping a matrix. */
    virtual ErrorFlag GetMatrixName( Matrix &A , std::string &name );

    /** Setup solver to be the default solver. This is used in case of failure. */
    virtual ErrorFlag GetDefaultSolver( Solver &solver ); 

};

}
/** Template implentations need the header, so we include them here */
#include "UserInterface.tpp"

#endif
