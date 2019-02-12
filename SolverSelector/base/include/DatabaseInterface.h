#ifndef SSDATABASEINTERFACE
#define SSDATABASEINTERFACE
//There are three types of hashes. Solver hash is a hash of the linear solver.
// Matrix hash (M) is a hash of the Matrix 
// Unique hash (U) is a hash of the solver and matrix combination and should be unique in the database. 

// The update status tells the database what to do when 
// a row is alrady in the database, Ignore means leave it as is,
// replace means replace it, and AVERAGE means take the average. A database might 
// not support all of these though. 

#include "Solvers.h"

namespace SolverSelecter { 

  /** Status for adding new rows to database. IGNORE -- if a row exists do nothing, REPACE -- if row
   * exists replace it, AVERAGE -- if a row exists, update the average value. */
  enum class UpdateStatus { IGNORE, REPLACE, AVERAGE } ;

  /** Class to interact with the database, whatever it may be */
  class DatabaseInterface : public SSBase { 
  public:
    
    
    bool initialized = false; /** has the database been initialized */
    bool reclassify = false;  /** do we need to reclassify the solvers in the database */
    UpdateStatus status;  /** what is the status for new rows */

    /* Constructor */
    DatabaseInterface ();
    
    /* Descructor */
    virtual ~DatabaseInterface() {};

    /** Create a database if it doesn't exist. */ 
    virtual ErrorFlag Initialize() = 0;

    /** Add a new row to the database */ 
    virtual int AddToDatabase(Solver &solver , 
                              std::map<std::string, double> &features,
                              std::map<std::string,double> &measurements,
                              UpdateStatus status) = 0;
    
    /** Add a solver classification to the database. */
    virtual ErrorFlag 
    AddClassificationToDatabase( int &hash, 
                                 std::map<std::string, bool > &classification, 
                                 UpdateStatus status ) = 0;
    
    /** Get the classification for the solver hash u */
    virtual ErrorFlag 
    GetUniqueClassification( int u, 
                             std::map<std::string, bool > &classification ) = 0;
    
    /** Get the measurement map for the row in database with hash $hash */
    virtual ErrorFlag 
    GetUniqueMeasurement( int &hash, 
                          std::map<std::string, double> &mvalues ) = 0;
    
    /** Get the feature set for the row in the database with hash $hash */
    virtual ErrorFlag 
    GetUniqueFeatures( int &hash, 
                       std::map<std::string, double> &features) = 0;

    /** Get the solver in the database with hash $hash */
    virtual ErrorFlag
    GetUniqueSolver( int &hash , 
                     Solver &solver ) =0;    

    /** Get the hash of the solver for a given row ($hash) in the database. */
    virtual ErrorFlag 
    GetUniqueSolverInRow( int &hash , 
                          int &solverhash ) =0;  

    /** Get a list of the unique solver hashes in the database */
    virtual ErrorFlag 
    GetUniqueSolverList( std::vector< int > &solvers ) = 0; 

    /** Get a list of the rows in the main database ( i.e, the measurements ) */
    virtual ErrorFlag 
    GetUniqueHashList( std::vector< int > &solvers ) = 0;
    
    /** Get a list of the features in the database. */
    virtual ErrorFlag 
    GetFeatureLabels( std::vector< std::string > &labels ) = 0;
    
    /** Get a list of the classifications in the database. */
    virtual ErrorFlag 
    GetClassificationLabels( std::vector< std::string > &clas ) = 0; 

    /** Get the map of matrices to solver hashes from the database. That is, 
     * get a list of solver hashes for each matrix. */
    virtual ErrorFlag 
    GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap, 
                          int matrix=-1) =0 ; 
    
    /** Pack up the database for input into a machine learning model. */
    virtual ErrorFlag 
    GetMachineLearningData(std::vector <std::string>&CNames,std::vector< int > &row_ids , 
                           std::vector< int > &solvers_labels, 
                           std::vector< std::string > &features_labels ,
                           std::vector< std::string > &classification_labels,
                           std::vector< int > &solvers_data,  
                           std::vector< std::vector< double > > &feature_data,
                           std::vector< std::vector< bool > > &classification_data );   

    /** Pre desctruction destructor ( finish up all MPI calls here */ 
    virtual ErrorFlag Finalize();
  
    /** Add a row to the database */
    virtual ErrorFlag 
    AddRow( Solver &solver, 
            std::map< std::string, double > &features, 
            std::map<std::string, double> &measurements);
    
    /* Add a classification to database */
    virtual ErrorFlag 
    AddClassification( int  &urow , 
                       std::map< std::string, double > &cvalues,
                       std::string CName, UpdateStatus status ) ; 
    
    /* Classify all the solvers in the database */
    virtual ErrorFlag 
    ClassifySolvers(  std::string MName, double MValue, std::string CName  ); 
    
};

}

#endif




