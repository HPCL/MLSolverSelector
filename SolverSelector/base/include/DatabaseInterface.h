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

  enum class UpdateStatus { IGNORE, REPLACE, AVERAGE } ;

  class DatabaseInterface : public SSBase { 
  public:
    
    bool initialized = false; 
    bool reclassify = false;  
    UpdateStatus status; 

    DatabaseInterface ();
    
    virtual ~DatabaseInterface() {};

    //Create a database if it doesn't exist. 
    virtual ErrorFlag  
    Initialize() = 0;

    virtual int
    AddSolverToDatabase( Solver &solver , 
                         UpdateStatus status) = 0;
    
    virtual int 
    AddMatrixToDatabase( std::map< std::string, double > &features, 
                         UpdateStatus status ) = 0;
    
    virtual int  
    AddRowToDatabase( int &solver, 
                      int &marix, 
                      std::map<std::string, double> &measurements, 
                      UpdateStatus status) = 0;

    virtual ErrorFlag 
    AddClassificationToDatabase( int &hash, 
                                 std::map<std::string, bool > &classification, 
                                 UpdateStatus status ) = 0;
    
    virtual ErrorFlag 
    GetUniqueClassification( int u, 
                             std::map<std::string, bool > &classification ) = 0;
    
    virtual ErrorFlag 
    GetUniqueMeasurement( int &hash, 
                          std::map<std::string, double> &mvalues ) = 0;
    virtual ErrorFlag 
    GetUniqueFeatures( int &hash, 
                       std::map<std::string, double> &features) = 0;

    virtual ErrorFlag
    GetUniqueSolver( int &hash , 
                     Solver &solver ) =0;    

    virtual ErrorFlag 
    GetUniqueSolverInRow( int &hash , 
                          int &solverhash ) =0;  

    virtual ErrorFlag 
    GetUniqueSolverList( std::vector< int > &solvers ) = 0; 

    virtual ErrorFlag 
    GetUniqueHashList( std::vector< int > &solvers ) = 0;
    
    virtual ErrorFlag 
    GetFeatureLabels( std::vector< std::string > &labels ) = 0;
    
    virtual ErrorFlag 
    GetClassificationLabels( std::vector< std::string > &clas ) = 0; 

    virtual ErrorFlag 
    GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap, 
                          int matrix=-1) =0 ; 

    virtual ErrorFlag 
    Finalize();

    virtual ErrorFlag 
    AddRow( Solver &solver, 
            std::map< std::string, double > &features, 
            std::map<std::string, double> &measurements);
    
    virtual ErrorFlag 
    AddClassification( int  &urow , 
                       std::map< std::string, double > &cvalues,
                       UpdateStatus status ) ; 
    
    virtual ErrorFlag 
    ClassifySolvers( std::map< std::string , double > &bvalues, 
                     int just_matrix_in_row=-1 ); 
    
    virtual ErrorFlag 
    GetMachineLearningData(std::vector< int > &row_ids , 
                           std::vector< int > &solvers_labels, 
                           std::vector< std::string > &features_labels ,
                           std::vector< std::string > &classification_labels,
                           std::vector< int > &solvers_data,  
                           std::vector< std::vector< double > > &feature_data,
                           std::vector< std::vector< bool > > &classification_data );
  };   
 
}

#endif




