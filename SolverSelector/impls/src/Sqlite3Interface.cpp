
#ifdef WITH_SQLITE3 
#include "Sqlite3Interface.h"

namespace SolverSelecter
{

/**
    * Constructor
    **/
SqlDatabase::SqlDatabase() 
{
    cm = " , ";
    lb =" ( ";
    rb =" ) ";
    qq ="\"";
    sc =";";
    matrix_table = "Matrix";
    meas_table = "Measurements";
    solver_table = "Solvers";
    clas_table = "Classification";

    AddParameter("Database" , {} , "Database Name" );
     
}

SqlDatabase::~SqlDatabase() {
  if ( db != NULL ) 
    sqlite3_close( db );
}

int
SqlDatabase::GetHash( std::string hash_str )
{
    int hash;
    HashString( hash_str.c_str(), hash );
    return hash;
}

int
SqlDatabase::SQLExecute( std::string sql )
{
    if ( ! sql.empty() )
    {
        sql_result.clear();
        sql_column.clear();
        rc = sqlite3_exec(db,sql.data(), callback, this, &zErrMsg);

        if ( rc != SQLITE_OK )
        {
            fprintf(stderr, " %s \n " , sql.c_str() );
            fprintf(stderr, " SQL ERROR: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
    }
    return error_flag;
}

int
SqlDatabase::callback(void *ptr /**< no idea */,
                                 int argc  /**< no idea */,
                                 char *argv[] /**< no idea */,
                                 char *azColName[] /**< no idea */)
{

    SqlDatabase * thisptr = reinterpret_cast<SqlDatabase* >(ptr);

    std::vector<std::string> result, column;
    thisptr->sql_result.push_back({});
    thisptr->sql_column.push_back({});
    for ( int i = 0; i < argc ; i++ )
    {
        thisptr->sql_result.back().push_back( ( argv[i] == NULL ) ? std::string() : argv[i] );
        thisptr->sql_column.back().push_back( ( azColName[i] == NULL ) ? std::string() : azColName[i] );
    }
    return 0;
}

ErrorFlag
SqlDatabase::Initialize() 
{
    std::string database_name = GetParameter("Database");

    std::cout << " Opening the database " << database_name << std::endl; 

    rc = sqlite3_open( database_name.c_str(), &db );
    if( rc )
      throw SSException("Can't open database at: " + database_name );
    
    /* Create SQL table if it doesn't exist  */
    std::ostringstream meas, mat, sol, clas ;
    meas << "CREATE TABLE IF NOT EXISTS "  << meas_table << "("
         << " UID INTEGER PRIMARY KEY NOT NULL,"
         << " SOLVER INT NOT NULL,"
         << " MATRIX INT NOT NULL"
         << " );";
    SQLExecute( meas.str() );

    // Matrix table -- features added as needed
    mat << "CREATE TABLE IF NOT EXISTS " << matrix_table << "("
        << " UID INTEGER PRIMARY KEY NOT NULL"
        << " );";
    SQLExecute( mat.str() );

    // Solvers table -- lookup based on hash
    sol << "CREATE TABLE IF NOT EXISTS " << solver_table << "("
        << " UID INTEGER PRIMARY KEY NOT NULL,"
        << " NAME TEXT NOT NULL,"
        << " PRECONDITIONER TEXT,"
        << " PSTRING TEXT"
        << " );";
    SQLExecute( sol.str() );

    //Classification table ( added as needed. Rows are measure table ids
    clas << "CREATE TABLE IF NOT EXISTS " << clas_table << "("
         << " UID INTEGER PRIMARY KEY NOT NULL "
         << ");";
    SQLExecute( clas.str() );
    
    initialized = true;
    return error_flag;
}


std::vector<std::pair<std::string, std::string> > 
SqlDatabase::GetColumnInfo(std::string table)
{

    if ( columns_updated.find(table) == columns_updated.end() )
    {
        std::vector<std::pair<std::string, std::string > > v;
        columns_updated[table] = std::make_pair( true, v);
    }

    if ( columns_updated[table].first == false )
        return columns_updated[table].second;

    columns_updated[table].first = false;
    std::vector< std::pair< std::string, std::string >  > column_types;

    std::ostringstream ss;
    ss << "PRAGMA table_info(" << table << ");";
    SQLExecute(ss.str());

    for ( unsigned int j = 0; j < sql_result.size(); j++)
    {
        std::string name, type;
        for ( unsigned int  i = 0; i < sql_result[j].size(); i++ )
        {

            if ( sql_column[j][i] == "name")
                name = sql_result[j][i];
            else if ( sql_column[j][i] == "type" )
                type = sql_result[j][i];
        }
        column_types.push_back( std::make_pair( name, type ) );
    }

    columns_updated[table] = std::make_pair( false, column_types );

    return column_types;
}

void
SqlDatabase::SQL_AddColumn( std::string table, std::string name, std::string type ,std::string def )
{
    std::ostringstream c;
    std::vector< std::pair<std::string,std::string> > cnames = GetColumnInfo(table);

    if ( ! CheckColumnExists( cnames, name ) )
    {
        c << "ALTER TABLE " << table << " ADD COLUMN " << name << " " << type ;

        if ( ! def.empty() )
            c << " DEFAULT " << def  ;
        c << sc;
        SQLExecute( c.str() );
        columns_updated[table].first = true;
    }
    return;
}

bool
SqlDatabase::CheckColumnExists( std::vector<std::pair<std::string, std::string> > &columns, std::string cname )
{

    auto it = std::find_if( columns.begin(), columns.end() ,
                            [&](const std::pair< std::string, std::string> &element)
    {
        return element.first == cname ;
    } );
    return it != columns.end();
}

bool
SqlDatabase::CheckRowExists( int hash, std::string table )
{
    std::ostringstream oss;
    oss << "SELECT UID FROM " << table << " WHERE UID=" << hash << sc;
    SQLExecute( oss.str() ) ;
    return sql_result.size() > 0 ;
}

std::vector<std::string>
SqlDatabase::AddNewColumns( std::vector< std::pair< std::string, int >  > &map, std::string table, std::string type )
{
    std::vector< std::pair<std::string,std::string> > cnames = GetColumnInfo(matrix_table);
    std::vector<std::string> new_cols;
    for ( auto &it : map )
    {
        if ( !CheckColumnExists( cnames, it.first ) )
        {
            new_cols.push_back(it.first);
            SQL_AddColumn( table, it.first, type );
            if ( it.second == 2 )
            {
                SQL_AddColumn( table, "av_" + it.first , "INTEGER", "0" );
            }
        }
    }
    return new_cols;
}

ErrorFlag
SqlDatabase::NewRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string,int >> &map )
{

    std::ostringstream pre,post;
    pre << "INSERT OR REPLACE INTO " << table << lb << "UID" ;
    post << " VALUES " << lb << hash ;

    for ( auto &it : map )
    {
        if ( it.second.second == 0 )
        {
            pre << cm << it.first ;
            post << cm << qq << it.second.first << qq ;
        }
        else if ( it.second.second == 1 )
        {
            pre << cm << it.first ;
            post << cm << it.second.first ;
        }
        else
        {
            pre << cm << it.first << cm << "av_" + it.first ;
            post << cm << it.second.first << cm << 1 ;
        }
    }

    pre << rb;
    post << rb << sc ;
    SQLExecute( pre.str() + post.str() );
    return 0;
}

ErrorFlag
SqlDatabase::ReplaceRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string, int >> &map )
{
    std::ostringstream pre,post;
    pre << "UPDATE " << table << " SET " ;


    int i = 0;
    for ( auto &it : map )
    {

        if ( i++ > 0 )
            pre << cm ;

        if ( it.second.second == 0 )
            pre << it.first << "=" << qq << it.second.first << qq;
        else
            pre << it.first << "=" << it.second.first ;

        if ( it.second.second == 2 )
            pre << cm << "av_" + it.first << "= 1 "; // Reset the average counter
    }

    pre << " WHERE UID=" << hash << sc;
    if ( i > 0 )
        SQLExecute( pre.str() );
    return 0;
}

ErrorFlag
SqlDatabase::AverageRowInDatabase( int hash, std::string table, std::map< std::string, std::pair<std::string, int > > &map )
{

    std::ostringstream pre,post;
    pre << "UPDATE " << table << " SET " ;
    int i = 0;
    for ( auto &it : map )
    {
        if ( i++ > 0 ) pre << cm ;

        if ( it.second.second == 0 )
            pre << it.first << "=" << qq << it.second.first << qq;
        else if ( it.second.second == 1 )
            pre << it.first << "=" << it.second.first ;
        else
        {
            std::ostringstream counter, floatcounter, current;

            pre << it.first << "=" << lb << it.first << " * " << "av_" + it.first << "+" << it.second.first << rb
                << " / " << "CAST" << lb << "av_" + it.first << " AS REAL " << rb ;
            pre << cm << "av_" << it.first << "=" << "av_" + it.first << " + 1 " ;
        }
        pre << " WHERE UID=" << hash << sc;
    }
    return 0;
}

ErrorFlag
SqlDatabase::NewRowInDatabase( int hash, std::string table,
                                   std::map< std::string, std::pair<std::string, int > > &map,
                                   std::vector<std::string> &new_rows  )
{

    //In this case, pick out the new rows and call a replace on just them
    std::map< std::string, std::pair< std::string, int > > new_map;
    for ( auto &it : new_rows )
        new_map.insert( std::make_pair( it , map[it] ) );
    if ( new_map.size() > 0 )
        ReplaceRowInDatabase( hash, table, new_map ) ;
    return 0;
}


ErrorFlag
SqlDatabase::SetRowInDatabase( int hash, std::string table, std::string type,
                                   std::map<std::string, std::pair< std::string , int >> &map, UpdateStatus status )
{

    // Add any new columns to the database.
    std::vector< std::pair< std::string, int >  > cnames;
    for ( auto &it : map ) cnames.push_back( std::make_pair( it.first, it.second.second ) ) ;
    std::vector< std::string > new_cols = AddNewColumns( cnames, table, type);

    if ( !CheckRowExists( hash, table ) )
    {
        NewRowInDatabase( hash, table, map );
    }
    else if ( status == UpdateStatus::REPLACE )
    {
        ReplaceRowInDatabase( hash, table, map );
    }
    else if ( status == UpdateStatus::AVERAGE )
    {
        AverageRowInDatabase( hash, table, map );
    }
    else
    {
        NewRowInDatabase( hash, table, map, new_cols ) ;
    }
    return 0;
}

std::vector < std::pair<std::string , std::string> >
SqlDatabase::GetUniqueRow( int hash, std::string table )
{

    std::vector< std::pair<std::string, std::string> > map;
    std::ostringstream get;
    get << "SELECT * FROM " << table << " WHERE UID=" << hash << sc;
    SQLExecute( get.str() );
    if ( sql_result.size() > 0 )
    {
        for ( int i = 0; i < sql_result[0].size(); i++ )
            map.push_back( std::make_pair( sql_column[0][i] , sql_result[0][i] ) );
    }
    else 
    {
      std::string e = "Row " + std::to_string(hash) + "does not exist in table " + table ; 
      throw SSException(e); 
    }
    return map;
}

int
SqlDatabase::AddSolverToDatabase( Solver &solver , UpdateStatus status )
{

    std::set< std::string > parameters;
    std::string solvername, precondname, parameterstring;
    solver.GetSolverInfo(solvername,precondname,parameters);
    solver.GetSolverString(parameterstring);

    std::map< std::string , std::pair< std::string, int >> map;
    map["NAME"] = std::make_pair( solvername, 0 );
    map["PRECONDITIONER"] = std::make_pair( precondname, 0 );
    map["PSTRING"] = std::make_pair( parameterstring, 0 );
    int solver_hash = GetHash(parameterstring);
    SetRowInDatabase( solver_hash, solver_table, "TEXT", map, status );

    return solver_hash;
}

int
SqlDatabase::AddMatrixToDatabase( std::map< std::string, double > &features, UpdateStatus status )
{

    std::map<std::string, std::pair< std::string, int > > map;

    std::stringstream oss;
    for ( auto &it : features )
    {
        oss << it.first << it.second ;
        map[it.first] = std::make_pair( std::to_string(it.second) , 2 );
    }

    int matrix_hash = GetHash( oss.str() );

    SetRowInDatabase( matrix_hash, matrix_table, "REAL", map, status );
    return matrix_hash;
}

int 
SqlDatabase::AddRowToDatabase( int &rowhash, int &solver, int &matrix, std::map<std::string, double> &measurements, UpdateStatus status )
{


    std::map<std::string, std::pair< std::string, int > > map;
    map["SOLVER"] = std::make_pair( std::to_string(solver) , 1 );
    map["MATRIX"] = std::make_pair( std::to_string(matrix) , 1 );


    for ( auto &it : measurements )
    {
        map[it.first] = std::make_pair( std::to_string(it.second) , 2 );
    }

    SetRowInDatabase( rowhash, meas_table, "REAL", map, status );
    return rowhash;
}

int SqlDatabase::AddToDatabase( Solver &solver, std::map<std::string, double> &features, std::map<std::string,double> &measurements, UpdateStatus status) {

    int shash = AddSolverToDatabase( solver, status );
    int mhash = AddMatrixToDatabase( features, status); 
    
    std::string rowHashStr;
    solver.GetSolverString(rowHashStr);
    std::stringstream oss;
    oss << rowHashStr ; 
    for ( auto &it : features )
        oss << it.first << it.second ;
    int rowHash = GetHash(oss.str() );
    
    
    return AddRowToDatabase( rowHash, shash, mhash, measurements, status);

}



ErrorFlag 
SqlDatabase::AddClassificationToDatabase( int &hash, std::map<std::string, bool > &classification, UpdateStatus status )
{

    std::map<std::string, std::pair<std::string, int > > map;
    for ( auto &it : classification )
    {
        map[it.first] = std::make_pair( std::to_string( it.second ) , 1 ) ;
    }
    SetRowInDatabase( hash, clas_table, "INTEGER", map, status ) ;
    return 0;
}


ErrorFlag 
SqlDatabase::GetUniqueClassification( int u, std::map<std::string, bool > &classification )
{

    std::vector<std::pair<std::string,std::string>> map = GetUniqueRow( u, clas_table );
    // Want 1:n because 0 is the uid.
    for ( int i = 1; i < map.size(); i++ )
        classification.insert( std::make_pair( map[i].first , std::atoi(map[i].second.c_str()) ) );
    return 0;
}

ErrorFlag
SqlDatabase::GetUniqueMeasurement( int &hash, std::map<std::string, double> &mvalues )
{
    std::vector<std::pair<std::string,std::string>> map = GetUniqueRow( hash, meas_table );
    // uid, solver, matrix, measurements
    for ( int i = 3; i < map.size() ; i++ )
        mvalues.insert( std::make_pair( map[i].first , std::atof(map[i].second.c_str()) ) );
    return 0;
}

ErrorFlag 
SqlDatabase::GetUniqueFeatures( int &hash, std::map<std::string, double> &features)
{
    //Here hash is the hash for the row in the measurements table. So, first need to extract
    //the matrix hash

    std::string matrix_hash = GetElement( hash, meas_table, "MATRIX" );
    int m = std::atoi(matrix_hash.c_str());
    std::vector<std::pair<std::string,std::string>> map = GetUniqueRow(m , matrix_table );
    // uid, name, features
    for ( int i = 1; i < map.size() ; i++ )
    {
        features.insert( std::make_pair( map[i].first , std::atof(map[i].second.c_str()) ) );
    }
    return 0;
}

std::string
SqlDatabase::GetElement( const int &hash, std::string table, std::string element )
{
    std::ostringstream oss;
    oss << " SELECT " << element << " FROM " << table << " WHERE UID=" << hash << sc ;
    SQLExecute( oss.str() ) ;
    if ( sql_result.size() > 0 && sql_result[0].size() > 0 )
        return sql_result[0][0] ;
    else 
    {
        std::string e = std::to_string(hash) + "does not exist in table or " + element \
                        + " does not exist in row " ;
        throw SSException(e);
    }    
}

ErrorFlag 
SqlDatabase::GetUniqueSolverInRow( int &hash , int &solverhash )
{
   try {
      solverhash = std::atoi(GetElement( hash, meas_table, "SOLVER" ).c_str());
   } catch (SSException e) {
      std::string ee = std::to_string(hash) + " not found in " + meas_table ;  
      throw SSException(ee);
   }
   return 0;
}

ErrorFlag
SqlDatabase::GetUniqueSolver( int &hash, Solver &solver )
{
    try { 
        std::string solverstr = GetElement( hash, solver_table, "PSTRING" );
        solver.ParseSolverString( solverstr.c_str() );
    } catch (SSException e) {
        std::string ee = "Unable to build Solver: "; 
        throw SSException(ee);
    }
    return 0;
}

ErrorFlag
SqlDatabase::GetUniqueSolverList( std::vector< int > &solvers )
{
    std::ostringstream oss;
    oss << " SELECT UID FROM " << solver_table << sc;
    SQLExecute( oss.str() );
    for ( auto it : sql_result )
        solvers.push_back( std::atoi( it[0].c_str() ) );
    return 0;
}

ErrorFlag
SqlDatabase::GetUniqueHashList( std::vector< int > &solvers )
{
    std::ostringstream oss;
    oss << "SELECT UID FROM " << meas_table << sc;
    SQLExecute(oss.str());
    for ( auto &it : sql_result )
        solvers.push_back( std::atoi ( it[0].c_str() ) ) ;
    return 0;
}

// List all features and classifications available in the database.
ErrorFlag 
SqlDatabase::GetFeatureLabels( std::vector< std::string > &labels )
{
    std::vector<std::pair<std::string, std::string>> info = GetColumnInfo(matrix_table);
    for ( int i = 1; i < info.size(); i++ )
    {
        if ( info[i].first.substr(0,3).compare("av_") != 0 )
            labels.push_back( info[i].first );
    }
    return 0;
}

ErrorFlag 
SqlDatabase::GetClassificationLabels( std::vector< std::string > &clas )
{
    std::vector< std::pair< std::string, std::string > > info = GetColumnInfo(clas_table);
    for ( int i = 1; i < info.size(); i++ )
        clas.push_back( info[i].first );
    return 0;
}

ErrorFlag 
SqlDatabase::GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap, int matrix_in_row )
{

    int matrix = -1;
    if ( matrix_in_row >= 0 )
    {
        std::ostringstream os;
        os << "SELECT MATRIX FROM " << meas_table << " WHERE UID = " << matrix_in_row << sc;
        SQLExecute( os.str() );
        int matrix = std::atoi(sql_result[0][0].c_str());
        std::ostringstream oss;
        oss << "SELECT MATRIX, UID FROM " << meas_table << " WHERE MATRIX=" << matrix << sc ;
        SQLExecute( oss.str() );
    }
    else
    {
        std::ostringstream oss;
        oss << "SELECT MATRIX,UID FROM " << meas_table << sc;
        SQLExecute( oss.str() );
    }

    //Sqlresult is a list of matrices->solver in the database.
    for ( auto &it : sql_result )
    {
        auto itt = solvermap.find(std::atoi(it[0].c_str()));
        if ( itt == solvermap.end() )
            solvermap[ std::atoi(it[0].c_str()) ] = { std::atoi(it[1].c_str()) };
        else
            solvermap[std::atoi(it[0].c_str())].push_back( std::atoi(it[1].c_str() ) ) ;
    }
    return 0;
}
/*  This function was used to read Kanikas arff files into the existing sql database. It is 
 *  probably not needed anymore. 
int SqlDatabase::WriteDatabaseFromArff(std::string arffFile) {
 
  // Step 1. Get a map of solver hashes from kanika to SS solvers 
  std::ifstream file ("/home/boneill/Documents/RNET_Development/RNET/MLNEAMS/MLdata/solvers.txt");
  std::string value;
  int count = 0;
  int solver = 0;
  std::map<int, std::string > solverMap;

  std::vector<std::string> result;
  while (!file.eof() ) {
      
      result.clear();
      std::getline( file, value);
      
      StringSplit(value, ",", result); 
      if ( result.size() == 2 )
        solverMap[std::atoi(result[0].c_str())] = result[1];
  } 
  
  // Step 2. Get a map of index to Feature name. 
  std::map<int, std::string > featureMap;
  std::ifstream file1 ("/home/boneill/Documents/RNET_Development/RNET/MLNEAMS/MLdata/MOOSE/RS1/RS1_petsc_mfree_MOOSE_artemis_p1_30_30nnz>10000.arff");
  int index = 0;
  int class_index = -1;
  int solver_index = -1;
  while ( !file1.eof() ) {
    result.clear();
    value = "";
    std::getline( file1, value );
    result.clear();
    StringSplit(value, " ", result);
    if ( value.size() == 0 ) {

    }
    else if ( result[0].substr(0,5) == "@data")  {
       break;
    } else if ( result[1].substr(0,6) == "solver" ) {
        solver_index = index++; 
        
    } else if ( result[1].substr(0,5) == "class") {
        class_index = index++; 
    }
    else if ( result[0].substr(0,5) == "@attr" && result[1].substr(0,1) != "" ) {
        featureMap[index++] = result[1];
    }  

  }

  std::vector< Solver > solvers;
  std::vector< std::map< std::string , double >> features;
  std::vector< std::map< std::string, double >> measurements;
  std::vector< std::map< std::string, bool >> classifications;


  // This BEGIN TRANSACTION stuff should be implemented in the whole 
  // thing. Basically, SQL writes to disk once per transaction. Writing
  // to disk is very slow becaise SQL does it in very SAFE way. So, to 
  // get any kind of performance, we need to group our calls together into
  // one big transaction. 1000 per transaction is the max, so 500 seems like
  // a good compromize. SQL does about 1 transaction per second -- so, 145000
  // lines should take ~10 mins or something. 
  int progress = 0;
  int cprogress = 0;
  while ( !file1.eof() ) {
    result.clear();
    std::getline(file1, value);
    StringSplit(value, ",",result);
  
    if ( progress == 0 ) {
       std::string s = "BEGIN TRANSACTION"; 
       SQLExecute(s);
       progress++;
    } else if ( progress == 200 ) {
       std::cout << progress*(cprogress+1) << " records converted " ; 
       std::string s = "END TRANSACTION" ;
       SQLExecute(s);
       progress = 0;
       cprogress++;
    } else {
      progress++;
    }
    std::string solverstr = solverMap[ std::atoi(result[solver_index].c_str()) ];
    Solver solver( solverstr );

    std::map< std::string, double> feature_map;
    for ( auto f : featureMap ) {
       feature_map[f.second] = std::atof(result[f.first].c_str());  
    }
    
    // Fake a measurement map that will always give us the given classification 
    bool classification = ( result[class_index] == "bad" ) ? false : true ;     
    std::map< std::string, double > measurementMap_unavailable;
    measurementMap_unavailable["CPUTime"] = ( classification ) ? 0.003 : 1e300 ; 
    int rowHash = AddToDatabase( solver, feature_map, measurementMap_unavailable, UpdateStatus::REPLACE );
 
    std::map<std::string, bool> classMap;
    classMap["CPUTime"] = classification;
    int classHash = AddClassificationToDatabase( rowHash, classMap, UpdateStatus::REPLACE );
    
  }
  
  if ( progress != 0 ) {
       std::string s = "END TRANSACTION" ;
       SQLExecute(s);
   } 
	 return 0;

}*/

}
#endif

