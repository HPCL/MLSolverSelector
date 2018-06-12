
#ifndef _SQLLITE3A_SS_HEADER
#define _SQLLITE3A_SS_HEADER


/* To use this, make sure you edit the makefile_inc to point to the sqlite3 install directory. */
#include "typedefs.h"
#include "DatabaseInterface.h"
#include "sqlite3.h"

/**
 * This is the SQLITE3 implimentation of the database base class.
 **/

class _SS_SqlDatabase : public _SS_DatabaseInterface
{

public:

    sqlite3 *db;                                                         /**< pointer to the sqlite3 database. */
    char *zErrMsg = 0;                                                   /**< stores the latest sqlite3 error message*/
    int rc;                                                              /**< stores the latest sqlite3 error code */
    std::vector < std::vector< std::string > > sql_result;               /**< stores the latest sqlite3 result */
    std::vector < std::vector< std::string > > sql_column;               /**< stores the latest sqlite3 result */

    std::string matrix_table, meas_table, clas_table, solver_table;       /**< the name of the table in the database */
    std::vector< std::string > column_names;                             /**< list of all the column names. */
    std::vector< std::string > column_types;                             /**< list of the the column types  */
    std::string cm,lb,rb,qq,sc ;                                             /**< some convienece strings used throughout  */

    /**
     * Constructor
     **/
    _SS_SqlDatabase( const std::string &_database_name /**< name of the database */) : _SS_DatabaseInterface(_database_name) {
        cm = " , ";
        lb =" ( ";
        rb =" ) ";
        qq ="\"";
        sc =";";
        matrix_table = "Matrix";
        meas_table = "Measurements";
        solver_table = "Solvers";
        clas_table = "Classification";
        std::cout << database_name << " is the databae \n" ; 
    }
    
    int GetHash( std::string hash_str )
    {
        int hash;
        _SS_Utils::HashString( hash_str.c_str(), hash );
        return hash;
    }
   
    int SQLExecute( std::string sql ) 
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
            return _SS_error_flag;
     }

     static int callback(void *ptr /**< no idea */,
                         int argc  /**< no idea */,
                         char **argv /**< no idea */,
                         char **azColName /**< no idea */) {
        
        _SS_SqlDatabase * thisptr = reinterpret_cast<_SS_SqlDatabase* >(ptr);
        
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


    virtual int Initialize() override
    {
        rc = sqlite3_open( database_name.c_str(), &db );
        if( rc )
            fprintf(stderr, "Can't open database at %s : %s\n", database_name.c_str(), sqlite3_errmsg(db));

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
            << " UID INTEGER PRIMARY KEY NOT NULL,"
            << " NAME TEXT NOT NULL"
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

        return _SS_error_flag;
    }
    
    std::vector<std::pair<std::string, std::string> > GetColumnInfo(std::string table) {
        
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
        return column_types;
    }
  
    std::ostringstream SQL_AddColumn( std::string table, std::string name, std::string type ,std::string def="" ) {
        std::ostringstream c;
        std::vector< std::pair<std::string,std::string> > cnames = GetColumnInfo(table);       
       
        if ( ! CheckColumnExists( cnames, name ) )
        {
          c << "ALTER TABLE " << table << " ADD COLUMN " << name << " " << type ;

          if ( ! def.empty() )
             c << " DEFAULT " << def  ;
          c << sc; 
          return c;
        }
    }

    bool CheckColumnExists( std::vector<std::pair<std::string, std::string> > &columns, std::string cname ) {
        
       auto it = std::find_if( columns.begin(), columns.end() , 
            [&](const std::pair< std::string, std::string> &element){ return element.first == cname ; } );
       return it != columns.end(); 
    }

    bool CheckRowExists( int hash, std::string table ) {
      std::ostringstream oss;
      oss << "SELECT UID FROM " << table << " WHERE UID=" << hash << sc;
      SQLExecute( oss.str() ) ;
      return sql_result.size() > 0 ; 
    }

    std::vector<std::string> AddNewColumns( std::vector< std::pair< std::string, int >  > &map, std::string table, std::string type )
    {
        std::vector< std::pair<std::string,std::string> > cnames = GetColumnInfo(matrix_table);
        std::vector<std::string> new_cols; 
        for ( auto &it : map ) {
           if ( !CheckColumnExists( cnames, it.first ) ) {
              new_cols.push_back(it.first);
              std::ostringstream oss = SQL_AddColumn( table, it.first, type );
              SQLExecute( oss.str() ); 
              if ( it.second == 2 ) { 
                std::ostringstream oss1 = SQL_AddColumn( table, "av_" + it.first , "INTEGER", "0" );
                SQLExecute( oss1.str() );
              }
           }
        }
        return new_cols;
   }

   int NewRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string,int >> &map ) {
        
        std::ostringstream pre,post; 
        pre << "INSERT OR REPLACE INTO " << table << lb << "UID" ;
        post << " VALUES " << lb << hash ; 
        
        for ( auto &it : map ) {
           if ( it.second.second == 0 ) { 
             pre << cm << it.first ; 
             post << cm << qq << it.second.first << qq ; 
           }
           else if ( it.second.second == 1 ) { 
             pre << cm << it.first ; 
             post << cm << it.second.first ;
           }
           else { 
             pre << cm << it.first << cm << "av_" + it.first ;
             post << cm << it.second.first << cm << 1 ; 
           }        
        }
        
        pre << rb; 
        post << rb << sc ; 
        SQLExecute( pre.str() + post.str() );
   }
 
   int ReplaceRowInDatabase( int hash, std::string table, std::map<std::string, std::pair<std::string, int >> &map ) {
      std::ostringstream pre,post;
      pre << "UPDATE " << table << " SET " ;
     
      
      int i = 0;
      for ( auto &it : map ) {
         
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

    int AverageRowInDatabase( int hash, std::string table, std::map< std::string, std::pair<std::string, int > > &map ) {

      std::ostringstream pre,post;
      pre << "UPDATE " << table << " SET " ;
      int i = 0; 
      for ( auto &it : map ) { 
        if ( i++ > 0 ) pre << cm ; 

        if ( it.second.second == 0 ) 
          pre << it.first << "=" << qq << it.second.first << qq; 
        else if ( it.second.second == 1 ) 
          pre << it.first << "=" << it.second.first ;
        else {
            std::ostringstream counter, floatcounter, current; 
            
            pre << it.first << "=" << lb << it.first << " * " << "av_" + it.first << "+" << it.second.first << rb 
                << " / " << "CAST" << lb << "av_" + it.first << " AS REAL " << rb ;
            pre << cm << "av_" << it.first << "=" << "av_" + it.first << " + 1 " ;  
        } 
        pre << " WHERE UID=" << hash << sc;
      }
      return 0;
    }

    int NewRowInDatabase( int hash, std::string table, 
                          std::map< std::string, std::pair<std::string, int > > &map,
                          std::vector<std::string> &new_rows  ) {

      //In this case, pick out the new rows and call a replace on just them 
      std::map< std::string, std::pair< std::string, int > > new_map;
      for ( auto &it : new_rows ) 
          new_map.insert( std::make_pair( it , map[it] ) );
      if ( new_map.size() > 0 )
        ReplaceRowInDatabase( hash, table, new_map ) ; 
      return 0;
    }
   

   int SetRowInDatabase( int hash, std::string table, std::string type,
                         std::map<std::string, std::pair< std::string , int >> &map, UpdateStatus status ) {

      // Add any new columns to the database. 
     std::vector< std::pair< std::string, int >  > cnames; 
     for ( auto &it : map ) cnames.push_back( std::make_pair( it.first, it.second.second ) ) ; 
     std::vector< std::string > new_cols = AddNewColumns( cnames, table, type);

      if ( !CheckRowExists( hash, table ) ) {  
          NewRowInDatabase( hash, table, map );
      } else if ( status == UpdateStatus::REPLACE ) {
         ReplaceRowInDatabase( hash, table, map );
      } else if ( status == UpdateStatus::AVERAGE ) { 
          AverageRowInDatabase( hash, table, map );
      } else { 
          NewRowInDatabase( hash, table, map, new_cols ) ;  
      }
      return 0;
   }

   std::vector < std::pair<std::string , std::string> > 
   GetUniqueRow( int hash, std::string table ) {
   
     std::vector< std::pair<std::string, std::string> > map;
     std::ostringstream get;
     get << "SELECT * FROM " << table << " WHERE UID=" << hash << sc;
     SQLExecute( get.str() ); 
     
      // Sql result should have one line in it;
     for ( int i = 0; i < sql_result[0].size(); i++ ) 
         map.push_back( std::make_pair( sql_column[0][i] , sql_result[0][i] ) );
     return map;
  }
      
  int AddSolverToDatabase( _SS_Solver &solver , UpdateStatus status ) {

      std::set< std::string > parameters;
      std::string solvername, precondname, parameterstring;
      solver.GetSolverInfo(solvername,precondname,parameters);
      solver.GetSolverString(parameterstring);
      
      std::map< std::string , std::pair< std::string, int >> map;
      map["NAME"] = std::make_pair( solvername, 0 );
      map["PRECONDITIONER"] = std::make_pair( precondname, 0 );
      map["PSTRING"] = std::make_pair( parameterstring, 0 );
      
      SetRowInDatabase( GetHash(parameterstring), solver_table, "TEXT", map, status );

      return 0;
   }

   int AddMatrixToDatabase( std::string matrix, std::map< std::string, double > &features, UpdateStatus status ) {

      std::map<std::string, std::pair< std::string, int > > map;
      map["NAME"] = std::make_pair( matrix, 0 );
      for ( auto &it : features ) {
        map[it.first] = std::make_pair( std::to_string(it.second) , 2 ); 
      }
      SetRowInDatabase( GetHash(matrix), matrix_table, "REAL", map, status );
      return 0;
   } 

   int AddRowToDatabase( _SS_Solver &solver, std::string &matrix, std::map<std::string, double> &measurements, UpdateStatus status ) {

      std::string solver_string, solvername;
      solver.GetSolverString(solver_string);

      std::map<std::string, std::pair< std::string, int > > map;
      map["SOLVER"] = std::make_pair( std::to_string(GetHash(solver_string)) , 1 );
      map["MATRIX"] = std::make_pair( std::to_string(GetHash(matrix)) , 1 );
     
      
      for ( auto &it : measurements ) {
        map[it.first] = std::make_pair( std::to_string(it.second) , 2 ); 
      }
      SetRowInDatabase( GetHash(solver_string + matrix), meas_table, "REAL", map, status );
      return 0;
   } 
      

  int AddClassificationToDatabase( int &hash, std::map<std::string, bool > &classification, UpdateStatus status ) {
      
    std::map<std::string, std::pair<std::string, int > > map;
    for ( auto &it : classification ) {
       map[it.first] = std::make_pair( std::to_string( it.second ) , 1 ) ; 
    }
    SetRowInDatabase( hash, clas_table, "INTEGER", map, status ) ; 
    return 0;
  }
  


  int GetUniqueClassification( int u, std::map<std::string, bool > &classification ) {
    std::vector<std::pair<std::string,std::string>> map = GetUniqueRow( u, clas_table );
    // Want 1:n because 0 is the uid. 
    for ( int i = 1; i < map.size(); i++ ) 
      classification.insert( std::make_pair( map[i].first , std::atoi(map[i].second.c_str()) ) );
    return 0;
  }

  int GetUniqueMeasurement( int &hash, std::map<std::string, double> &mvalues ) { 
    std::vector<std::pair<std::string,std::string>> map = GetUniqueRow( hash, meas_table );
    // uid, solver, matrix, measurements  
    for ( int i = 3; i < map.size() ; i++ ) 
      mvalues.insert( std::make_pair( map[i].first , std::atof(map[i].second.c_str()) ) );
    return 0;
  }
  
  int GetUniqueFeatures( int &hash, std::map<std::string, double> &features) {
      std::vector<std::pair<std::string,std::string>> map = GetUniqueRow( hash, matrix_table );
      // uid, name, features   
      for ( int i = 2; i < map.size() ; i++ ) 
         features.insert( std::make_pair( map[i].first , std::atof(map[i].second.c_str()) ) );
    return 0;
  }
  
  std::string GetElement( const int &hash, std::string table, std::string element ) {
    std::ostringstream oss;
    oss << " SELECT " << element << " FROM " << table << " WHERE UID=" << hash << sc ;
    SQLExecute( oss.str() ) ;
    return sql_result[0][0] ; 
  }

  int GetUniqueMatrix( int &hash , std::string &matrix ) {
      //Get the matrix associated with this row 
    std::string matrix_hash = GetElement( hash, meas_table, "MATRIX" );
    matrix = GetElement( std::atoi( matrix_hash.c_str() ) , matrix_table, "NAME" ) ;
    return 0; 
  }
 
  int GetUniqueSolver( int &hash , _SS_Solver &solver ) {
    std::string solverhash = GetElement( hash, meas_table, "SOLVER" );
    std::string solverstr = GetElement( std::atoi(solverhash.c_str()), solver_table, "PSTRING" ); 
    solver.ParseSolverString( solverstr.c_str() );
    return 0;
  }    
  
  int GetUniqueHashList( std::vector< int > &solvers ) { 
     std::ostringstream oss;
     oss << "SELECT UID FROM " << meas_table << sc; 
     for ( auto &it : sql_result ) 
        solvers.push_back( std::atoi ( it[0].c_str() ) ) ; 
     return 0;
  }

  // List all features and classifications available in the database. 
  int GetFeatureLabels( std::vector< std::string > &labels ) { 
    std::vector<std::pair<std::string, std::string>> info = GetColumnInfo(matrix_table);
    for ( int i = 2; i < info.size(); i++ ) 
       labels.push_back( info[i].first );
  }

  int GetClassificationLabels( std::vector< std::string > &clas ) {
    std::vector< std::pair< std::string, std::string > > info = GetColumnInfo(clas_table);
    for ( int i = 1; i < info.size(); i++ )
       clas.push_back( info[i].first );
    return 0;
  } 

  int GetMatrixToUniqueMap( std::map< int, std::vector<int> > &solvermap ) {
    std::ostringstream oss;
    oss << "SELECT MATRIX,UID FROM " << meas_table ;
    SQLExecute( oss.str() );
     
    //Sqlresult is a list of matrices->solver in the database. 
    for ( auto &it : sql_result ) {
       auto itt = solvermap.find(std::atoi(it[0].c_str()));
       if ( itt == solvermap.end() )  
            solvermap[ std::atoi(it[0].c_str()) ] = { std::atoi(it[1].c_str()) };
       else 
            solvermap[std::atoi(it[0].c_str())].push_back( std::atoi(it[1].c_str() ) ) ; 
    }
    return 0;
  }
} ;

#endif
