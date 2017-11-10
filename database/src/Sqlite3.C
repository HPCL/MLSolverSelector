
#ifdef WITH_SQLITE3
/** \file Sqlite3.C 
 * \brief Implimentations for database integration with sqlite3. This file is included at the end
 * of the Sqlite3.h file. That allows us to remain header only, but still seperate implimentations
 * and header files.  
 **/
#include "Sqlite3.h"

using namespace Database;

_SS_DataBaseSql::_SS_DataBaseSql( const std::string &_database_name ) 
  : _SS_DataBaseBase( _database_name)
{
  cm = " , ";
  lb =" ( ";
  rb =" ) ";
  qq ="\"";   
  data_table = "DATA";

}

_SS_ErrorFlag
_SS_DataBaseSql::Initialize()
{
  rc = sqlite3_open( database_name.c_str(), &db );  
  if( rc ) 
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));

  /* Create SQL table if it doesn't exist  */
  std::ostringstream ss;
  ss << "CREATE TABLE IF NOT EXISTS "  << data_table << "(" 
     << " _S_UID INT PRIMARY KEY NOT NULL,"
     << " _S_HASH INT NOT NULL,"
     << " _S_NAME TEXT NOT NULL," 
     << " _S_MATRIX TEXT NOT NULL,"
     << " _S_PRECONDITIONER TEXT,"
     << " _S_PARAMETERS TEXT"   
     << " );";
  //create the table //
  SQLExecute( ss.str() );

  //Get the table infomation 
  GetColumnNames();
  return _SS_error_flag;
}

_SS_ErrorFlag
_SS_DataBaseSql::Finalize()
{
   sqlite3_close(db); 
   return _SS_error_flag;
}

_SS_ErrorFlag _SS_DataBaseSql::GetColumnNames( )
{
  std::ostringstream ss;
  ss << "PRAGMA table_info(" << data_table << ");";
  SQLExecute(ss.str());
  column_names.clear(); 
  for ( unsigned int j = 0; j < sql_result.size(); j++)
  {
    for ( unsigned int  i = 0; i < sql_result[j].size(); i++ )
    {
      if ( sql_column[j][i] == "name")
        column_names.push_back( sql_result[j][i] );
      if ( sql_column[j][i] == "type")
        column_types.push_back( sql_result[j][i] );
    }
  }
  return _SS_error_flag;
  
}

_SS_ErrorFlag _SS_DataBaseSql::GetSolver( int &hash, _SS_Solver &solver )
{

  std::stringstream ss;
  ss << " SELECT _S_PARAMETERS FROM " << data_table << " WHERE _S_HASH=" << hash; 
  SQLExecute(ss.str());
  std::string sstring;
  if ( sql_result.size() >0  && sql_result[0].size() > 0 )
  {
    /* There is probably more than one row, but we only need one  */
    sstring = sql_result[0][0];
    solver.ParseSolverString( sstring ); 
  }
  else
    std::cout << " Error getting the solver : hash does not exist \n ";

  return _SS_error_flag;
}

_SS_ErrorFlag _SS_DataBaseSql::GetSolverHashList( std::vector< int > & shashes )
{
  std::stringstream ss;
  ss << " SELECT DISTINCT _S_HASH FROM " << data_table << ";";
  SQLExecute(ss.str());
  
  shashes.clear();
  for ( auto it : sql_result )
  {
      shashes.push_back( std::stoi(it[0]) );
  }
  
  return _SS_error_flag;
}


_SS_ErrorFlag _SS_DataBaseSql::CheckColumnExists( const std::string &column_name , bool &exists)
{
    exists = std::find(column_names.begin(), column_names.end(), column_name) != column_names.end();       
    return _SS_error_flag;
}


_SS_ErrorFlag _SS_DataBaseSql::AddColumn( const std::string &prefix, const std::string &name, const std::string &type  )
{
    bool already_exists;
    CheckColumnExists( prefix + name, already_exists );

    if ( ! already_exists )
    {
      std::ostringstream ss;
      ss << "ALTER TABLE "<< data_table << " ADD COLUMN " << prefix + name << " " << type << ";";
      SQLExecute(ss.str());      
      GetColumnNames();
    }

    return _SS_error_flag; 
}

_SS_ErrorFlag _SS_DataBaseSql::CheckRowExists( const int &uhash, int &exists )
{
  /* Function to check if a row is already in the database. Each row has a hash
   * of the string "matrix solver precondtioner pname pvale ptype pname pvalue ..... */
  
  std::ostringstream oss ;
  oss << "SELECT * FROM " << data_table << " WHERE _S_UID=" << uhash <<";";
  SQLExecute(oss.str());

  if ( sql_result.size() > 0 )
    exists = 1;
  else
    exists = 0;

  return _SS_error_flag;
}

_SS_ErrorFlag _SS_DataBaseSql::ModifyRow( const int &uhash, 
                                          const _SS_Measurements &measurements,
                                          const _SS_Features &features)
{
 
  /* Modify the values in the row. Currently, this is a straight replace. Could 
   * also do a "average" or whatever else */
  bool exists;
  int run = 0;
 
  std::map< std::string, double> mparams, fparams;
  std::string fprefix, mprefix, ftype;
  

  fprefix = "_F_";
  mprefix = "_M_";
  measurements.Get(mparams );
  features.Get(fparams );

  std::ostringstream oss;
  oss << "UPDATE " << data_table << " SET ";
  
  if (mparams.size() > 0 )
  {
    for ( auto it : mparams )
    {

      CheckColumnExists( mprefix + it.first  , exists);
      if (exists)
      {
        oss << mprefix << it.first << "=" << it.second << cm;
        run++;
      } 
    }
  }
  if ( fparams.size() > 0 )
  {
    for ( auto it : fparams )
    {
      CheckColumnExists( fprefix + it.first  , exists);
      if (!exists)
          AddColumn(fprefix, it.first, "REAL");
      oss << fprefix << it.first << "=" << it.second << cm;
      run++;
       
    }
  }
   
  if ( run > 0 )
  {
    long pos = oss.tellp();  
    oss.seekp(pos-3);
    oss << " WHERE _S_UID="<< uhash << ";";
    SQLExecute(oss.str());
  }
  return _SS_error_flag;
} 

_SS_ErrorFlag _SS_DataBaseSql::AddRow( const _SS_Solver &solver, 
                                       const std::string &matrix,
                                       const _SS_Features &features,
                                       const _SS_Measurements &measurement )
{
   int uhash, shash,  exists;
   GetHash( solver, matrix, uhash, shash );
   CheckRowExists( uhash, exists );
   if ( exists )
   {
     ModifyRow( uhash, measurement, features);
   }
   else
   {  
     std::set< std::string > parameters;
     std::string solvername, precondname, parameterstring;
     solver.GetSolverInfo(solvername,precondname,parameters);
     solver.GetSolverString(parameterstring);

     /* Now that that is all sorted, we can add the row */
     std::ostringstream ss_pre, ss_names, ss_values, ss_post;
     ss_pre << "INSERT INTO " << data_table << " ";
     
     ss_names  << lb << "_S_UID"            << cm
                     << "_S_HASH"           << cm
                     << "_S_NAME"           << cm 
                     << "_S_PRECONDITIONER" << cm 
                     << "_S_MATRIX"         << cm 
                     << "_S_PARAMETERS"     << cm;
     
     ss_values << " VALUES ";
     ss_values << lb << qq << uhash << qq << cm
                     << qq << shash << qq << cm
                     << qq << solvername << qq << cm 
                     << qq << precondname << qq << cm
                     << qq << matrix << qq << cm
                     << qq << parameterstring << qq << cm;

     std::map < std::string, double > feats,meas; 
     features.Get( feats );
     std::string fprefix = "_F_";
     for ( auto key : feats )
     {
       AddColumn(fprefix, key.first, "REAL");
       ss_names    << fprefix << key.first << cm ;
       ss_values << key.second << cm ;      
     }
      
     /* Loop over all the measurements */
     bool exists1;
     std::string mprefix = "_M_";
     measurement.Get( meas );
     for ( auto &key : meas )
     {
       AddColumn(mprefix, key.first, "REAL" );
       ss_names    << mprefix << key.first << cm ;
       ss_values << key.second << cm ;             
       CheckColumnExists( "_C_" + key.first  , exists1);
       if ( exists1 )
       {
          ss_names << "_C_" << key.first << cm;
          ss_values << 0 << cm ; 
       }
     }
    
     CheckColumnExists( "_R_restriction", exists1); 
     if (exists1)
     {
        ss_names << "_R_restriction" << " ) ";
        ss_values << 0 << " ) " ;
     }
     else
     {   
       /* Delete the last comma */
       long pos = ss_values.tellp();  ss_values.seekp(pos-3);
       long pos1 = ss_names.tellp();  ss_names.seekp(pos1-3);
       ss_values << "  )   "; 
       ss_names << "  )   "; 
     }
     std::string sql = ss_pre.str() + ss_names.str() + ss_values.str() + ";" ;
     
     SQLExecute( sql );
     
   }
  return _SS_error_flag;
}

_SS_ErrorFlag
_SS_DataBaseSql::RestrictSolvers( std::vector < std::pair < std::string, std::string > > &sset )
{

  /* Add a column to the data base */
  AddColumn( "_R_", "restriction", "INT" );
  if ( sset.size() > 0 )
  {
    /* Set everything to zero */
    std::stringstream oss;
    oss << " UPDATE " << data_table << " SET _R_restriction=0; ";
    SQLExecute(oss.str());

    /* Set every row that matches to one */
    for ( auto it : sset ) 
    {
      std::stringstream sss;
      sss << " UPDATE " << data_table << " SET _R_restriction = 1 WHERE _S_Solver=" << it.first << "\"" 
          << " AND _S_PRECONDITIONER=" << it.second << "\";";
      SQLExecute(sss.str());
    }
  }
  else
  {
    std::stringstream oss;
    oss << " UPDATE " << data_table << " SET _R_restriction=1 ; ";
    SQLExecute(oss.str());
  }

  return _SS_error_flag;
}
  
  _SS_ErrorFlag 
_SS_DataBaseSql::ImportData( std::vector< std::pair< std::string , std::string > > & sset, 
                             std::vector< std::vector < std::string > > &data , 
                             std::vector<std::string> &ret_column_names,
                             std::vector<int> &feature_or_label )
{
  
  /* restrict the solvers that can be used */  
  RestrictSolvers(sset);

  /* Set the attributes and build sql string*/
  std::stringstream ss;
  
  data.clear();
  ss << "SELECT _S_HASH" ;  
  
  ret_column_names.clear();
  ret_column_names.push_back("HASH");
  feature_or_label.push_back(0);

  /* Go through once and get the features */
  for ( auto c : column_names )
  {
    std::string prefix( c.substr(0,3) );
    if ( prefix == "_F_" )
    {
       ss << cm << c ;
       ret_column_names.push_back(c.substr(3));
       feature_or_label.push_back(0);
    }
  }  
   
  /* Go through another time and get the labels. */ 
  for ( auto c : column_names )
  {
    std::string prefix( c.substr(0,3) );
    if ( prefix == "_C_" )
    {
       ss << cm << c ;
       ret_column_names.push_back(c.substr(3));
       feature_or_label.push_back(1);
    }
  }  
  

  ss << " FROM " << data_table << " WHERE _R_restriction = 1 ";
  SQLExecute(ss.str());
  data = sql_result; 


  return _SS_error_flag;
}


/* This is the sql callback function. It takes the call back and sets the 
 * sql_result and sql_column vectors with it */
int _SS_DataBaseSql::callback(void *ptr, int argc, char **argv, char **azColName)
{
  _SS_DataBaseSql * thisptr = reinterpret_cast<_SS_DataBaseSql* >(ptr);
  std::vector< std::string > result, column;

  for ( int i=0; i< argc; i++)
  {
    char buff[10000], buff1[10000];
    snprintf(buff,sizeof(buff), "%s", argv[i] );
    snprintf(buff1,sizeof(buff1), "%s", azColName[i] );
    result.push_back(buff);
    column.push_back(buff1);
  }
  thisptr->sql_result.push_back( result );    
  thisptr->sql_column.push_back( column );
  
  return 0;
}

_SS_ErrorFlag
_SS_DataBaseSql::SQLExecute( const std::string &sql )
{ 
  sql_result.clear();
  sql_column.clear();
  rc = sqlite3_exec(db,sql.data(), callback, this, &zErrMsg);
  
  if ( rc != SQLITE_OK )
  {
    fprintf(stderr, " SQL ERROR: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }
  return _SS_error_flag;
}

_SS_DataBaseSql::~_SS_DataBaseSql( )
{
}

/* ***** Need a better way to do this. Also, not sure if classification process is correct or not.  */
_SS_ErrorFlag
_SS_DataBaseSql::ClassifySolvers( const _SS_Measurements &measure )
{
  /* Classify the solvers based on the percentage from minimum */
  /* Step 1. Get the minimum value for each Matrix */
  std::map< std::string, double > keys;
  bool exists;
  measure.Get(keys);

  if ( keys.size() > 0 )
  {
    std::ostringstream oss;
    oss << "SELECT _S_MATRIX " << cm ;
    for ( auto key : keys )
    {
      CheckColumnExists("_M_" + key.first , exists );
      if (exists)
        oss << "min" << lb << "_M_" << key.first << rb << cm ;   
    }
    long pos1 = oss.tellp();  oss.seekp(pos1-3);
    oss << " FROM " << data_table << " GROUP BY _S_MATRIX ;";       
    SQLExecute(oss.str());
    // Return is < < mat, min0,min1,min2... > < mat1,min0,min1,min2...> ... > 

    std::vector< std::vector< std::string > > sql_copy = sql_result;
    std::string matrix;
    double mvalue;
    std::map<std::string, double> badvals;
    measure.GetBad(badvals);

    for ( auto mvals : sql_copy )
    {  
        matrix = mvals[0];
        std::ostringstream oss; 
        oss << " UPDATE " << data_table << " SET " ;
        int i = 0;
        for ( auto key : keys )
        {
          
          AddColumn("_C_", key.first, "INT" );
          mvalue = badvals[key.first];
          i++;
          double value = ( mvals[i] == "(null)" ) ? 0 : std::stof(mvals[i])*(1+mvalue) ;
          oss << "_C_" << key.first << " =( CASE WHEN _M_" << key.first << " < " << value << " THEN 1 ELSE 0 END )";
          oss << cm ; 
        }  
        long pos1 = oss.tellp();  oss.seekp(pos1-3);
        oss << " WHERE _S_MATRIX=" << qq << matrix << qq << ";"; 
        SQLExecute(oss.str());
    }
  }
    return _SS_error_flag;
}    

#endif 
