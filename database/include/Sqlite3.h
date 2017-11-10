#ifndef _SQLLITE3_SS_HEADER
#define _SQLLITE3_SS_HEADER

#ifdef WITH_SQLITE3

/* To use this, make sure you edit the makefile_inc to point to the sqlite3 install directory. */
#include "DataBase.h"
#include "sqlite3.h"

namespace Database
{
/**
 * This is the SQLITE3 implimentation of the database base class. 
 **/

class _SS_DataBaseSql : public _SS_DataBaseBase
{

public:
  
 sqlite3 *db;                                                         /**< pointer to the sqlite3 database. */
 char *zErrMsg = 0;                                                   /**< stores the latest sqlite3 error message*/
 int rc;                                                              /**< stores the latest sqlite3 error code */
 std::vector < std::vector< std::string > > sql_result;               /**< stores the latest sqlite3 result */
 std::vector < std::vector< std::string > > sql_column;               /**< stores the latest sqlite3 result */

 std::string data_table;                                              /**< the name of the table in the database */ 
 std::vector< std::string > column_names;                             /**< list of all the column names. */
 std::vector< std::string > column_types;                             /**< list of the the column types  */
 std::string cm,lb,rb,qq;                                             /**< some convienece strings used throughout  */
                                                  
  
  /**
   * Constructor
   **/
  _SS_DataBaseSql( const std::string &_database_name /**< name of the database */); 
  
  /**
   *  Destuctor 
   **/
  virtual ~_SS_DataBaseSql( );
  
  /************ Base Class implimentations ******************************************************/


  /**
   * Initialize the database. In this case, this connects to the sqlite3 database, and creates 
   * a new table, if the table does not already exist. The column_names variable is also updated here.
   **/
  virtual _SS_ErrorFlag Initialize() override ; 
 
  /**
   * FInalize: Closes the connection to the sqlite3 database. 
   **/
  virtual _SS_ErrorFlag Finalize() override ; 
  
  /**
   * Classify the solvers 
   **/
  virtual _SS_ErrorFlag ClassifySolvers( const _SS_Measurements &measure) override ;
  
  /**
   *  Import the data (override)
   **/
  virtual _SS_ErrorFlag ImportData( std::vector< std::pair< std::string , std::string > > & sset, 
                             std::vector< std::vector < std::string > > &data , 
                             std::vector<std::string> &column_names,
                             std::vector<int> &feature_or_label ) override;

  /**
   *  Add a new row (see base class for details)
   **/
  virtual _SS_ErrorFlag AddRow( const _SS_Solver &solver, 
                                const std::string &matrix,
                                const _SS_Features &features,
                                const _SS_Measurements &measurement ) override;

  /**
   * Get the solver for the hash (override, see base class for details) 
   **/
  virtual _SS_ErrorFlag GetSolver( int &shash , _SS_Solver &solver ) override ;
  
  /**
   * Return a list of all solvers in the database (see base class for details)
   **/
  virtual _SS_ErrorFlag GetSolverHashList( std::vector< int > &solver_hashes ) override ;

  /**
   * Sqite3 keeps an up to date database, so no write to file is neccesary (override) 
   */
  virtual _SS_ErrorFlag WriteToFile() override {return _SS_error_flag; } 

  /************ Sub class functions ******************************************************/
  
  /**
   *  Executute the sql statment
   **/
  _SS_ErrorFlag SQLExecute(const std::string &sql /**< string representing a sql query */ );
  
  /**
   *  Default callback function for sqlite3 result. This takes the sqlite3 result and stores it in 
   *  the sql_result and sql_column variables 
   **/
  static int callback(void *ptr /**< no idea */,
                      int argc  /**< no idea */, 
                      char **argv /**< no idea */, 
                      char **azColName /**< no idea */);


  /**
   * Return the column names 
   **/
  virtual _SS_ErrorFlag GetColumnNames(); 
  /**
   *  Check if the column exists 
   **/
  virtual _SS_ErrorFlag CheckColumnExists( const std::string &column_name /**< name of the column to find */, 
                                           bool &exists /**< output, true if column exists */ );
  /**
   *  Add a new column to the database. This implimentation uses a prefix system to determine between 
   *  the parameters, features, measurements, classifications, etc. 
   **/
  virtual _SS_ErrorFlag AddColumn( const std::string &prefix /**< prefix for the column */, 
                                   const std::string &feature /**< name of the column */, 
                                   const std::string &type /**< type of data stored in this column */ );  
  /**
   *  Check if a row exists in the database. 
   **/
  virtual _SS_ErrorFlag CheckRowExists( const int &uhash /**< uhash to search for */, 
                                        int &exists /**< output, true if the row exists */ );
  /**
   *  Modify a row if it already exists. This adds the row if it doesn't exist, and replaces the row 
   *  if it does. Might add support for averaging rows in the future if it proves useful. 
   **/
  virtual _SS_ErrorFlag ModifyRow(const int &uhash /**< hash for the row in question */, 
                                  const _SS_Measurements &param /**< updated measurements */, 
                                  const _SS_Features &features /**< updated features */);
  /**
   * Restrict the solvers that can be used in the database. This creates a new column in the database called 
   * _R_restriction with a initial value of 1. If a row has a solver precond pair that is not in the 
   * list, then that value is set 0, and the data from that row is not imported. This is called during the 
   * importData function. 
   **/
  virtual _SS_ErrorFlag RestrictSolvers( std::vector < std::pair<  std::string, std::string > > &sset /**< list of allowed solver/precond */) ;


};
}
#endif
#endif 
