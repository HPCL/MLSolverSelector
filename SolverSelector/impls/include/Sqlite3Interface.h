#ifndef SQLITE3INTERFACE_H
#define SQLITE3INTERFACE_H

#ifdef WITH_SQLITE3 

#include "DatabaseInterface.h"
#include "sqlite3.h"

namespace SolverSelecter { 

class SqlDatabase : public DatabaseInterface
{

public:

    std::vector< std::string > column_names; 
    std::vector< std::string > column_types; 
    std::string cm,lb,rb,qq,sc ;          
    std::string matrix_table, meas_table, clas_table, solver_table;  
    std::vector < std::vector< std::string > > sql_result;           
    std::vector < std::vector< std::string > > sql_column;           
    std::map<std::string, std::pair< bool, std::vector<std::pair<std::string,std::string> > > > columns_updated;
    
    sqlite3 *db = NULL;                                                    
    char *zErrMsg = 0;                                               
    int rc;                                                          

    SqlDatabase();
    
    virtual ~SqlDatabase(); 


    int GetHash( std::string hash_str );
      
    int GetRowHash( Solver &solver, std::map<std::string, double> matrix );

    static int callback(void *ptr, int argc, char *argv[], char *azColName[] ) ;
    
    int SQLExecute( std::string sql );

    virtual ErrorFlag Initialize() override;

    std::vector<std::pair<std::string, std::string> > GetColumnInfo(std::string table);

    void SQL_AddColumn( std::string table, std::string name, std::string type ,std::string def="" );

    bool CheckColumnExists( std::vector<std::pair<std::string, std::string> > &columns, std::string cname );


    bool CheckRowExists( int hash, std::string table );

    std::vector<std::string> AddNewColumns( std::vector< std::pair< std::string, int >  > &map, std::string table, std::string type );

    int AddToDatabase(Solver &solver, std::map<std::string, double> &features,  std::map<std::string,double> &measurements, UpdateStatus status) override;


    ErrorFlag NewRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string,int >> &map );

    ErrorFlag ReplaceRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string, int >> &map );

    ErrorFlag AverageRowInDatabase( int hash, std::string table, std::map< std::string, std::pair<std::string, int > > &map );

    ErrorFlag NewRowInDatabase( int hash, std::string table,
                          std::map< std::string, std::pair<std::string, int > > &map,
                          std::vector<std::string> &new_rows  );


    ErrorFlag SetRowInDatabase( int hash, std::string table, std::string type,
                          std::map<std::string, std::pair< std::string , int >> &map, UpdateStatus status );

    std::vector < std::pair<std::string , std::string> > GetUniqueRow( int hash, std::string table ) ;

    int AddSolverToDatabase( Solver &solver , UpdateStatus status );

    int AddMatrixToDatabase( std::map< std::string, double > &features, UpdateStatus status );

    int AddRowToDatabase( int &rowhash, int &solverhash, int &matrixHash, std::map<std::string, double> &measurements, UpdateStatus status );

    ErrorFlag AddClassificationToDatabase( int &hash, std::map<std::string, bool > &classification, UpdateStatus status );

    ErrorFlag GetUniqueClassification( int u, std::map<std::string, bool > &classification );

    ErrorFlag GetUniqueMeasurement( int &hash, std::map<std::string, double> &mvalues );

    ErrorFlag GetUniqueFeatures( int &hash, std::map<std::string, double> &features);

    std::string GetElement( const int &hash, std::string table, std::string element ) ;

    ErrorFlag GetUniqueSolverInRow( int &hash , int &solverhash );

    ErrorFlag GetUniqueSolver( int &hash, Solver &solver );

    ErrorFlag GetUniqueSolverList( std::vector< int > &solvers );

    ErrorFlag GetUniqueHashList( std::vector< int > &solvers );

    ErrorFlag GetFeatureLabels( std::vector< std::string > &labels ) ;

    ErrorFlag GetClassificationLabels( std::vector< std::string > &clas );

    ErrorFlag GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap, int matrix_in_row=-1 );
} ;

}
#endif

#endif
