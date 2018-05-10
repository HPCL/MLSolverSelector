#ifndef _SQLLITE3_SS_HEADER
#define _SQLLITE3_SS_HEADER


/* To use this, make sure you edit the makefile_inc to point to the sqlite3 install directory. */
#include "DataBase.h"
#include "sqlite3.h"

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
    _SS_DataBaseSql( const std::string &_database_name /**< name of the database */) : _SS_DataBaseBase(_database_name) {
        cm = " , ";
        lb =" ( ";
        rb =" ) ";
        qq ="\"";
        data_table = "DATA";

    }


    /**
     *  Destuctor
     **/
    virtual ~_SS_DataBaseSql( ) {};

    /************ Base Class implimentations ******************************************************/


    /**
     * Initialize the database. In this case, this connects to the sqlite3 database, and creates
     * a new table, if the table does not already exist. The column_names variable is also updated here.
     **/
    virtual _SS_ErrorFlag Initialize() override
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

      
    _SS_ErrorFlag SQLExecute( std::string sql ) 
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
    
    /**
     * FInalize: Closes the connection to the sqlite3 database.
     **/
    virtual _SS_ErrorFlag Finalize() override {
        sqlite3_close(db);
        return _SS_error_flag;
    }


    /**
     * Classify the solvers
     **/
    virtual _SS_ErrorFlag ClassifySolvers( const std::map<std::string, double> &mvalues ) override {
        /* Classify the solvers based on the percentage from minimum */
        /* Step 1. Get the minimum value for each Matrix */
        bool exists;

        if (mvalues.size() > 0 )
        {
            std::ostringstream oss;
            oss << "SELECT _S_MATRIX " << cm ;
            for ( auto key : mvalues )
            {
                CheckColumnExists("_M_" + key.first , exists );
                if (exists)
                    oss << "min" << lb << "_M_" << key.first << rb << cm ;
            }
            long pos1 = oss.tellp();
            oss.seekp(pos1-3);
            oss << " FROM " << data_table << " GROUP BY _S_MATRIX ;";
            SQLExecute(oss.str());
            // Return is < < mat, min0,min1,min2... > < mat1,min0,min1,min2...> ... >

            std::vector< std::vector< std::string > > sql_copy = sql_result;
            std::string matrix;

            for ( auto mvals : sql_copy )
            {
                matrix = mvals[0];
                std::ostringstream oss;
                oss << " UPDATE " << data_table << " SET " ;
                int i = 0;
                for ( auto key : mvalues )
                {

                    AddColumn("_C_", key.first, "INT" );
                    i++;
                    double value = ( mvals[i] == "(null)" ) ? 0 : std::stof(mvals[i])*(1+key.second) ;
                    oss << "_C_" << key.first << " =( CASE WHEN _M_" << key.first << " < " << value << " THEN 1 ELSE 0 END )";
                    oss << cm ;
                }
                long pos1 = oss.tellp();
                oss.seekp(pos1-3);
                oss << " WHERE _S_MATRIX=" << qq << matrix << qq << ";";
                SQLExecute(oss.str());
            }
        }
        return _SS_error_flag;
    }


    /**
     *  Import the data (override)
     **/
    virtual _SS_ErrorFlag ImportData( std::vector< std::pair< std::string , std::string > > & sset,
                                      std::vector< std::vector < std::string > > &data ,
                                      std::vector<std::string> &ret_column_names,
                                      std::vector<int> &feature_or_label ) override {


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


        ss << " FROM " << data_table;
        SQLExecute(ss.str());
        data = sql_result;

        return _SS_error_flag;
    }


    /**
     *  Add a new row (see base class for details)
     **/
    virtual _SS_ErrorFlag AddRow( const _SS_Solver &solver,
                                  const std::string &matrix,
                                  const _SS_features_map &fparams,
                                  const _SS_measurements_map &mparams ) override {
        
        int uhash, shash,  exists;
        GetHash( solver, matrix, uhash, shash );
        CheckRowExists( uhash, exists );
        if ( exists )
        {
            printf("sdfsdf\n");
            ModifyRow( uhash, mparams, fparams);

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


            std::string fprefix = "_F_";
            for ( auto key : fparams )
            {
                AddColumn(fprefix, key.first, "REAL");
                ss_names    << fprefix << key.first << cm ;
                ss_values << key.second << cm ;
            }

            /* Loop over all the measurements */
            bool exists1;
            std::string mprefix = "_M_";
            for ( auto &key : mparams )
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
                long pos = ss_values.tellp();
                ss_values.seekp(pos-3);
                long pos1 = ss_names.tellp();
                ss_names.seekp(pos1-3);
                ss_values << "  )   ";
                ss_names << "  )   ";
            }
            std::string sql = ss_pre.str() + ss_names.str() + ss_values.str() + ";" ;
            
            printf("%s\n", sql.c_str());
            SQLExecute( sql );

        }
        return _SS_error_flag;
    }


    /**
     * Get the solver for the hash (override, see base class for details)
     **/
    virtual _SS_ErrorFlag GetSolver( int &hash , _SS_Solver &solver ) override {

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


    /**
     * Return a list of all solvers in the database (see base class for details)
     **/
    virtual _SS_ErrorFlag GetSolverHashList( std::vector< int > &shashes ) override {
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


    /**
     * Sqite3 keeps an up to date database, so no write to file is neccesary (override)
     */
    virtual _SS_ErrorFlag WriteToFile() override {
        return _SS_error_flag;
    }


    virtual _SS_ErrorFlag GetFeatureNames( std::vector< std::string > &feature_names ) override {
        GetColumnNames();
        for ( auto it : column_names )
        {
            std::string prefix( it.substr(0,3) );
            if ( prefix == "_F_" )
                feature_names.push_back( it.substr(3,it.size()));
        }
        return _SS_error_flag;
    }


    /**
     *  Default callback function for sqlite3 result. This takes the sqlite3 result and stores it in
     *  the sql_result and sql_column variables
     **/
    static int callback(void *ptr /**< no idea */,
                        int argc  /**< no idea */,
                        char **argv /**< no idea */,
                        char **azColName /**< no idea */) {
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



    /**
     * Return the column names
     **/
    virtual _SS_ErrorFlag GetColumnNames() {
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

    /**
     *  Check if the column exists
     **/
    virtual _SS_ErrorFlag CheckColumnExists( const std::string &column_name /**< name of the column to find */,
            bool &exists /**< output, true if column exists */ )
    {
        exists = std::find(column_names.begin(), column_names.end(), column_name) != column_names.end();
        return _SS_error_flag;
    }

    /**
     *  Add a new column to the database. This implimentation uses a prefix system to determine between
     *  the parameters, features, measurements, classifications, etc.
     **/
    virtual _SS_ErrorFlag AddColumn( const std::string &prefix /**< prefix for the column */,
                                     const std::string &name /**< name of the column */,
                                     const std::string &type /**< type of data stored in this column */ ) {
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

    /**
     *  Check if a row exists in the database.
     **/
    virtual _SS_ErrorFlag CheckRowExists( const int &uhash /**< uhash to search for */,
                                          int &exists /**< output, true if the row exists */ ) {
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

    /**
     *  Modify a row if it already exists. This adds the row if it doesn't exist, and replaces the row
     *  if it does. Might add support for averaging rows in the future if it proves useful.
     **/
    virtual _SS_ErrorFlag ModifyRow(const int &uhash /**< hash for the row in question */,
                                    const _SS_measurements_map &mparams /**< updated measurements */,
                                    const _SS_features_map &fparams /**< updated features */) {

        /* Modify the values in the row. Currently, this is a straight replace. Could
         * also do a "average" or whatever else */
        bool exists;
        int run = 0;

        std::string fprefix = "_F_";
        std::string mprefix = "_M_";

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

};

#endif
