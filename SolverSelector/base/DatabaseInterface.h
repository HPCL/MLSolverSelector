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

enum class UpdateStatus { IGNORE, REPLACE, AVERAGE } ;

class _SS_DatabaseInterface { 
  public:
  bool reclassify = false; // Set this flag true to reclassify the database on finalization. 
  UpdateStatus status; 
  std::string database_name;

  _SS_DatabaseInterface ( std::string _database_name )  {
      database_name = _database_name;
  }

  // Add a row and matrix 
  virtual int Initialize() = 0;
  virtual int Finalize() {};
  virtual int AddSolverToDatabase( _SS_Solver &solver , UpdateStatus status) = 0;
  virtual int AddMatrixToDatabase( std::string matrix, std::map< std::string, double > &features, UpdateStatus status ) = 0;
  virtual int AddRowToDatabase( _SS_Solver &solver, std::string &matrix, 
                                std::map<std::string, double> &measurements, UpdateStatus status) = 0;
  virtual int AddClassificationToDatabase( int &hash, std::map<std::string, bool > &classification, UpdateStatus status ) = 0;
  virtual int GetUniqueClassification( int u, std::map<std::string, bool > &classification ) = 0;
  virtual int GetUniqueMeasurement( int &hash, std::map<std::string, double> &mvalues ) = 0;
  virtual int GetUniqueFeatures( int &hash, std::map<std::string, double> &features) = 0;
  virtual int GetUniqueMatrix( int &hash , std::string &matrix ) =0;
  virtual int GetUniqueSolver( int &hash , _SS_Solver &solver ) =0; 
  virtual int GetUniqueHashList( std::vector< int > &solvers ) = 0;
  virtual int GetFeatureLabels( std::vector< std::string > &labels ) = 0;
  virtual int GetClassificationLabels( std::vector< std::string > &clas ) = 0; 
  virtual int GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap ) =0 ; 

  int AddRow( _SS_Solver &solver, 
              std::string &matrix, 
              std::map< std::string, double > &features, 
              std::map<std::string, double> &measurements) 
  {  
    AddSolverToDatabase( solver, status  ) ;
    AddMatrixToDatabase( matrix, features, status );
    AddRowToDatabase( solver, matrix, measurements, status);
    return 0;
  }
  
  int AddClassification( int  &urow , std::map< std::string, double > &cvalues, UpdateStatus status ) {
      std::map<std::string, double> mvalues;
     std::map<std::string, bool >  bvalues; 
     GetUniqueMeasurement( urow , mvalues );
  
     for ( auto &it : cvalues ) {
        auto m = mvalues.find(it.first); 
        if ( m == mvalues.end() || m->second < 0.0 ) {
           bvalues[it.first] = false;
        } else {
           bvalues[it.first] = ( m->second < it.second ) ? true : false ;
        }
     }
     AddClassificationToDatabase( urow, bvalues, UpdateStatus::REPLACE  );      
  }

  int GetSolver( int hash , _SS_Solver &solver ) {
    GetUniqueSolver( hash, solver ); 
    return 0;
  }

  int ClassifySolvers( std::map< std::string , double > &bvalues ) {
    
      std::cout << " Classifying the Solvers \n";

      std::map< int, std::vector< int > > matrixmap; 
      GetMatrixToUniqueMap( matrixmap ) ;
      
      for ( auto matrix : matrixmap ) {
          std::map<std::string, double> cvalues; 
          for ( auto solver : matrix.second ) { 
            std::map<std::string, double> mvalues;
            GetUniqueMeasurement( solver , mvalues ) ; 
            for ( auto & it : mvalues ) 
              std::cout << it.first << " " << it.second << "\n";

            for ( auto meas : bvalues ) {  
                auto m = mvalues.find(meas.first);
                auto c = cvalues.find(meas.first);
                if ( m != mvalues.end() && m->second > 0 ) {
                    double cutoff = ( 1.0 + meas.second )*m->second ;
                    if ( c == cvalues.end() )
                      cvalues[meas.first] = cutoff; 
                    else if ( cutoff < c->second ) { 
                         c->second = cutoff; 
                    }
                }
            }
          }
          for ( auto solver : matrix.second ) {
            for ( auto &it : cvalues )                
              AddClassification( solver, cvalues, UpdateStatus::REPLACE ) ; 
          }
      }
  }

  int GetMachineLearningData(std::vector< int > &row_ids , 
                             std::vector< std::string > &features_labels ,
                             std::vector< std::string > &classification_labels, 
                             std::vector< std::vector< double > > &feature_data,
                             std::vector< std::vector< bool > > &classification_data ) {


    std::vector< int > rows; 
    GetUniqueHashList( rows );  // List of every row in the database.
    GetFeatureLabels(features_labels);
    GetClassificationLabels(classification_labels);

    for ( auto it : rows ) {

        row_ids.push_back( it ); 

        std::map< std::string, double > fvalues ; 
        std::map<std::string, bool > cvalues;
        std::vector< double > features; 
        std::vector< bool > classi; 
        GetUniqueFeatures( it, fvalues );  
        GetUniqueClassification( it, cvalues );
        
        for ( auto itt : features_labels ) { 
           auto f = fvalues.find(itt);
           if ( f != fvalues.end() ) 
             features.push_back( f->second ) ;
           else
             features.push_back( 1e300 ); 
        }
        feature_data.push_back(features);

        for ( auto itt : classification_labels ) { 
           auto f = cvalues.find(itt);
           if ( f != cvalues.end() ) 
             classi.push_back( f->second ) ;
           else
             classi.push_back( false ); 
        }
        classification_data.push_back(classi);
     }
     return 0;
  }
};     
#endif

