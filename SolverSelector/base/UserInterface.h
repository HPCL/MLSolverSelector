#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Solvers.h"
#include "MLWaff.h"
#include "MLBinary.h"
#include "Sqlite3.h"

template<typename Matrix,typename Vector>
class _SS_UserInterface
{
public:
   
    std::string database_name = "default.db"; 
    bool matrix_dump;
    bool use_ml = true;
    bool build_inline = false; 
    std::string input_file = "petsc.input";
      
    /* This is a virtual base class that must be overwritten
     * by the user of the code. It holds all the user defined functions
     * required by the solver selecter. */
    _SS_UserInterface() {} ;

    /** User impliments this solve function. The options chosen by the solver
    * selecter can be obtained by calling solvers.GetChosenOptions(...)  */
    virtual _SS_ErrorFlag SolveSystem( Matrix &A, Vector &x, Vector &b, _SS_Solver &solver ) = 0 ;

    
    /** Extract the features from the matrix */
    virtual _SS_ErrorFlag ExtractFeatures( Matrix &A, std::map<std::string, double> &fmap ) = 0;
   
    /** Start Measureing the measurements */
    virtual _SS_ErrorFlag StartMeasurements( Matrix &A, Vector &x, Vector &b ) = 0;

    /** Stop Measurements */
    virtual _SS_ErrorFlag StopMeasurements( Matrix &A, Vector &x, Vector &b, 
                                            std::map<std::string, double> &mmap ) = 0;

    virtual _SS_ErrorFlag ClassificationMeasurements( std::map<std::string, double> &class_values ) 
    {
        std::cout << "Classification Measurements not implimented\n";
        return _SS_error_flag; 
    }

    virtual _SS_ErrorFlag GetDataBaseImplimentation(std::shared_ptr< _SS_DataBaseBase > &database)
    {
        database.reset( new _SS_DataBaseSql(database_name) );     
        return _SS_error_flag;
    } 
    
    virtual _SS_ErrorFlag GetMachineLearningModel(std::shared_ptr< _SS_MachineLearning > &machine)
    {
        machine.reset( new _ML_Waffles() ) ;
        return _SS_error_flag;
    }    

    /** (optional) dump matrix to file **/
    virtual _SS_ErrorFlag DumpMatrix( Matrix &A ) 
    {
        return _SS_error_flag;
    }

    /**
     *  This is an optional function only required if using the BuildDataBaseFromFile option.
     **/
    virtual _SS_ErrorFlag InitMatrix( std::string filename,
                                      std::unique_ptr<Matrix> &A )
    {
        std::cout << "Init Matrix Not Implimented\n";
        return _SS_error_flag;
    }

    /** 
     * Initialize a vector based on the Matrix in AA 
     * */
    virtual _SS_ErrorFlag InitVector( const Matrix &AA,
                                      std::unique_ptr<Vector> &x )
    {
        std::cout << "Init Vector Not Implimented\n";
        return _SS_error_flag;
    }
    
    /** 
     * Initialize a vector based on the Matrix in AA 
     * */
    virtual _SS_ErrorFlag CopyVector( const Vector &AA,
                                      std::unique_ptr<Vector> &x )
    {
        std::cout << "Copy Vector Not Implimented\n";
        return _SS_error_flag;
    }
    
    /** Set a vector based on the values TODO **/
    virtual _SS_ErrorFlag SetVector( std::unique_ptr<Vector> &cloned_x, std::string type_ ) 
    {
        /* type is either "ones" for a vector of ones, or "rand" for a random vector. */
        std::cout << "SetVector Not Implimented\n";
        return _SS_error_flag;
    }

    /**
     * Free the pointers needed during building of the database.
     **/
    virtual _SS_ErrorFlag FreeMatrix( std::unique_ptr<Matrix> &A ) 
    {
        std::cout << "Free Matrix Not Implimented\n";
        return _SS_error_flag;
    }

    /** Override if some actions must be completed before the vector is destroyed. I.e, for Petsc,
     * we must vall VecDestroy before the pointer can be destroyed. */
    virtual _SS_ErrorFlag FreeVector( std::unique_ptr<Vector> &x ) 
    {
        std::cout << "Free Vector Not Implimented\n";
        return _SS_error_flag;
    }

    /** When building a database inline ( i.e., inside the code ) we need to name the matrix somehow. By default,
     * the name is just inline_current_time */
    virtual _SS_ErrorFlag GetMatrixName( Matrix &A , std::string &name )
    {
        name = "default_name" + std::to_string( std::rand() );
        return _SS_error_flag;
    }

    /* This should be set up to return the default solver */ 
    virtual _SS_ErrorFlag GetDefaultSolver( _SS_Solver &solver ) 
    {
        std::cout << "GetDefaultSolver Not Implimented\n";
        return _SS_error_flag;
    }


};

#endif

