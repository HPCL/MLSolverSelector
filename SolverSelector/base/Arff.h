
#ifndef _ARFFDB_SS_HEADER
#define _ARFFDB_SS_HEADER


/** I Need to get around to testing this. It compiles at least, but thats as far as I have gotton. 

/* To use this, make sure you edit the makefile_inc to point to the sqlite3 install directory. */
#include "DataBase.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>

class _SS_DataBaseArff : public _SS_DataBaseBase
{

public:

    bool write_to_file = false; // This flag gets set if the data_base is modified 
    bool overwrite_database = true; // Set this flag to overrite existing or add new 

    std::string title = "Solver Selector Feature Set";
    std::string sources = "Created by the RNET Solver Selector";
    std::string relation = "feature_set";
    std::vector< std::pair< std::string, std::string > > attributes;
    std::map< int,  std::map< std::string, std::string> > data; // { Hash -- > { feature --> feature_value } } 

    /**
     * Constructor
     **/
    _SS_DataBaseArff( const std::string &_database_name /**< name of the database */) 
        : _SS_DataBaseBase(_database_name)
    {

    }

    /**
     *  Destuctor
     **/
    virtual ~_SS_DataBaseArff( ) {};

    /************ Base Class implimentations ******************************************************/
    
    _SS_ErrorFlag ImportArffDatabase()
    {
        std::fstream file;
        file.open(database_name.c_str(), std::ios::in | std::ios::out | std::ios::app | std::ios::binary );

        std::string line, l;
        std::vector< std::string > result; 
        int status = 0;
        while (!file.eof())
        {
            std::getline(file,line); /* get the next line in the file */
            result.clear();
            _SS_Utils::StringSplit( line, " ", result );
            

            std::cout << result[0] << " sdfsdf" << std::endl; 

            if ( result.size() == 0 ) continue;           
            else if ( (result[0].substr(0,1)) == "%" ) continue;
            else if ( result[0] == "@RELATION" )
            {
                if ( status != 0 )
                {
                   std::cout << "1MISFORMED FILE __ TODO ABORT \n";
                   return 0;
                }
                else 
                {
                   relation = result[1] ;
                   status = 1; 
                }
            }
            else if ( result[0] == "@ATTRIBUTE" ) 
            {
                if ( status != 1 || result.size() != 3 )
                {
                   std::cout << "2MISFORMED FILE __ TODO ABORT \n";
                   return 0;
                }
                else 
                {
                   attributes.push_back( std::make_pair( result[1], result[2] ) ) ;
                   status = 1; 
                }
           }
           else if ( result[0] == "@DATA" ) 
           {
                if ( status != 1 || result.size() != 1 )
                {
                   std::cout << "3MISFORMED FILE __ TODO ABORT \n";
                   return 0;
                }
                else 
                {
                   status = 2; 
                }
           }
           else
           {
              // This is a data row. So, we want to split on commas-- 
             std::vector< std::string > cresult;
             _SS_Utils::StringSplit( line, ",", cresult );
             if ( cresult.size() != attributes.size() )
             { 
                  std::cout << "4MISFORMED FILE __ TODO ABORT \n";
                  return 0;
             }
             else
             { 
                std::map< std::string, std::string > row_map;
                
                int uhash = std::stoi(cresult[0]);
                for ( int i = 1; i < attributes.size(); i++ ) 
                {
                    row_map[attributes[i].first] = cresult[i] ;
                }
                data[uhash] = row_map;
             }
          }
        }

        for (auto it : attributes ) std::cout << it.first << " " << it.second << std::endl; 

        return _SS_error_flag;
    }        
                
    /**
     * Initialize the database. In this case, this connects to the sqlite3 database, and creates
     * a new table, if the table does not already exist. The column_names variable is also updated here.
     **/
    virtual _SS_ErrorFlag Initialize() override
    {
        struct stat buffer;   
        if (stat (database_name.c_str(), &buffer) == 0)
        {
            // The file already exists, so read it! 
            ImportArffDatabase( );
            if ( ! overwrite_database ) 
            {
              while ( stat ( database_name.c_str(), &buffer) == 0 )
                  database_name += "_";
              std::cout << " New filename will be " << database_name << std::endl;
              write_to_file = true;
            } 
        }
        else
        {
            // The database doesn't exist, so set up the default file attirbutes 
            attributes.push_back( std::make_pair( "UHASH" , "STRING" ) );            //Hash for this particular row 
            attributes.push_back( std::make_pair( "SHASH" , "STRING" ) );            //Hash for this particular solver
            attributes.push_back( std::make_pair( "SOLVER" , "STRING" ) );            //Name for this solver
            attributes.push_back( std::make_pair( "MATRIX" , "STRING" ) );          //Name of the matrix 
            attributes.push_back( std::make_pair( "PRECONDITIONER" , "STRING" ) );  //Name of the preconditioner 
            attributes.push_back( std::make_pair( "PARAMETERS" , "STRING" ) );      //Parameter string for this solver. 
        }     
        

        return _SS_error_flag;
    }

          
    /**
     * FInalize: This is called before destruction -- so write the database to file.
     **/
    virtual _SS_ErrorFlag Finalize() override
    {
        if ( ! write_to_file ) return 0 ;
      
        std::ofstream ofs(database_name.c_str(), std::ofstream::out); 
        ofs.setf(std::ios_base::scientific);
        ofs << "%" << title << "\n";
        ofs << "%" << sources << "\n";
        ofs << "@RELATION " << relation << "\n";

        for ( auto &it : attributes )
        {
            ofs << "@ATTRIBUTE " << it.first << " " <<  it.second << "\n" ; 
        }
        ofs << "\n\n@DATA\n" ; 

        for ( auto &it : data )
        {
            ofs << it.first ; 

            for ( auto &itt : attributes )
            {
                if (itt.first == "UHASH") // already added -- 
                    continue; 

                auto ittt = it.second.find(itt.first) ; 
                if ( ittt == it.second.end() ) 
                  ofs << ",?"; 
                else if ( itt.second == "STRING" )
                  ofs << "," << Stringify(ittt->second); 
                else 
                  ofs << "," << ittt->second ;
            }
            ofs << "\n";
        }

        return _SS_error_flag;
    }


    /**
     * Classify the solvers
     **/
    virtual _SS_ErrorFlag ClassifySolvers( const std::map<std::string, double> &mvalues ) override {
      //This is very slow -- FIXME -- Figure out a faster way to do this.  
        write_to_file = true;
        
        //Get the minimum of the measurement for each of the matricies in the database. 
        std::map< std::string , std::map < std::string, double > > minimums ; 
        for ( auto &row : data ) 
        {
             std::string matrix = row.second.find("MATRIX")->second;
             if ( minimums.find(matrix) == minimums.end() ) 
             {    
                  minimums[matrix] = mvalues;
                  for ( auto &m : minimums[matrix] ) 
                        m.second = std::numeric_limits<double>::max() ;
             }
             
             // This matrix is in there already 
             for ( auto &m : mvalues ) 
             {
                 auto rowm = row.second.find(m.first); 
                 if ( rowm != row.second.end() ) std::cout << "a" <<  rowm->second << "b" << std::endl; 
                 if ( rowm != row.second.end() && std::stod(rowm->second) < minimums[matrix][m.first] ) 
                    minimums[matrix][m.first] = std::stod(rowm->second); 
             }        
        }

        // Now go through and classify those matricies as good or bad. The classification will
        // go in a column called C_name. First step is to check if this classification row exists 
        // allready.
        std::map< std::string, bool > exists;
        for ( auto &m : mvalues ) 
        {
            bool exists = false;
            std::string cname = "C_" + m.first ; 
            for ( auto &at : attributes )
            {
                 if ( at.first == cname )
                 {
                      exists = true;
                      break;
                 }
            }
            if ( !exists ) 
            {
                attributes.push_back( std::make_pair( cname, "{0,1}" ) );
            }
         }
        
        for ( auto &row : data )
        { 
            std::string matrix = row.second.find("MATRIX")->second;
            auto matrix_min = minimums[matrix]; 
            for ( auto &m : mvalues )
            {
                auto measured = row.second.find(m.first); 
                if ( measured != row.second.end() )
                {  
                    std::string cname = "C_" + m.first;
                    row.second[cname] = ( std::stod(measured->second) < (1.0 + m.second )*matrix_min[m.first] ) ? "1" : "0" ; 
                }
            }
        }

        return _SS_error_flag;
    }

    

    /**
     *  Import the data (override) sset doesn't work > 
     **/
    virtual _SS_ErrorFlag ImportData( std::vector< std::vector < std::string > > &ddata ,
                                      std::vector<std::string> &ret_column_names,
                                      std::vector<int> &feature_or_label ) override 
    {
          //TODO -- This is a bit silly. The Machine learning file takes this info and builds a arff files
          //from it to get the features. 
          //data is a vector of vectors containing the data -- This is the rows 
          std::vector< bool > destring; 
          ret_column_names.clear();
          ret_column_names.push_back("SHASH"); 
          feature_or_label.push_back(0);

          for ( auto &it : attributes )
          {
              if ( it.second == "NUMERIC" )
              { 
                  feature_or_label.push_back(0);
                  ret_column_names.push_back(it.first);
                  std::cout << " Adding feature " << it.first << std::endl;
              }
              else if ( it.second != "STRING" )
              { 
                  feature_or_label.push_back(1);
                  std::cout << " Adding label " << it.first << std::endl;
                  ret_column_names.push_back(it.first);
              }
              else
                  std::cout << " Adding Nothing " << it.first << std::endl;
          }

          int idx = 0;
          for ( auto &it : data )
          {
            std::vector< std::string > row_data ;
            for ( auto &at : ret_column_names )
            {
                auto f = it.second.find(at);
                if ( f == it.second.end() ) 
                    row_data.push_back("?");
                else
                {
                  row_data.push_back( DeStringify(f->second));
                }
            }
            ddata.push_back(row_data);
          } 

        return _SS_error_flag;
    }


    /**
     *  Add a new row (see base class for details)
     **/
    virtual _SS_ErrorFlag AddRow( const _SS_Solver &solver,
                                  const std::string &matrix,
                                  const _SS_features_map &fparams,
                                  const _SS_measurements_map &mparams ) override {
        

        // Add a row, so we need to write to file. 
        write_to_file = true;

        // At the moment we just overwrite a row if it already exists. 

        int uhash, shash,  exists;
        GetHash( solver, matrix, uhash, shash );
        
        std::set< std::string > parameters;
        std::string solvername, precondname, parameterstring;
        solver.GetSolverInfo(solvername,precondname,parameters);
        solver.GetSolverString(parameterstring);
       
        std::map< std::string, std::string > row_data; 

        row_data["SHASH"] = std::to_string(shash); 
        row_data["SOLVER"] = solvername;
        row_data["PRECONDITIONER"] = precondname; 
        row_data["MATRIX"] = matrix;
        row_data["PARAMETERS"] = parameterstring;
        
        std::vector< std::string > cnames = GetColumnNames();
        for ( auto &it : fparams )
        {
            if ( std::find(cnames.begin(), cnames.end(), it.first) == cnames.end() )
                attributes.push_back(std::make_pair(it.first, "NUMERIC")); 
            row_data[it.first] = std::to_string( it.second );
        }
        for ( auto &it : mparams )
        {
            //Classify measurements as STRING because they are technically extra data not needed.  
            if ( std::find(cnames.begin(), cnames.end(), it.first) == cnames.end() )
                attributes.push_back(std::make_pair(it.first, "STRING")); 
            
            std::stringstream out; 
            out << std::scientific << it.second;
            row_data[it.first] = out.str();
        }

        data[uhash] = row_data;

        return _SS_error_flag;
    }


    /**
     * Get the solver for the hash (override, see base class for details)
     **/
    virtual _SS_ErrorFlag GetSolver( int &hash , _SS_Solver &solver ) override 
    {
        //Need to search through the database and extract the solver for that hash. The 
        //hash is the id for a particualr solver. 
        std::string hashs = std::to_string(hash);  
        for ( auto &it : data )
        {
           if ( it.second["SHASH"] == hashs )
           { 
               // Found it -- > Fixme -_ Need to generate a seperate solver map so we 
               // can do this faster. 
             std::string pstring = (it.second)["PARAMETERS"];
               solver.ParseSolverString( pstring ); 
               return 0;
           }
        }
        return _SS_error_flag;
    }


    /**
     * Return a list of all solvers in the database (see base class for details)
     **/
    virtual _SS_ErrorFlag GetSolverHashList( std::vector< int > &shashes ) override {
        std::set< int > sset; 
        for ( auto &it : data )
        {
            sset.insert( std::stoi( it.second["SHASH"] ) );
        }
        
        std::copy(sset.begin(), sset.end(), std::back_inserter(shashes));
        return _SS_error_flag;
    }

    virtual _SS_ErrorFlag GetFeatureNames( std::vector< std::string > &feature_names ) override {
         
        for ( auto &it : attributes )
        {
            if ( it.second == "NUMERIC" )
               feature_names.push_back(it.first);
        }
        return _SS_error_flag;
    }

    std::string Stringify(std::string s) { 
      return "'" + s + "'" ;
    }
    std::string DeStringify(std::string s ) {
      if ( s.substr(0, 1) == "'" )
          return s.substr(1, s.size()-2);
      else
          return s;
      }


    /**
     * Return the column names
     **/
    std::vector< std::string> GetColumnNames()
    {
        std::vector< std::string > n;
        for ( auto &it : attributes )
           n.push_back( it.first );
        return n;
    }
      
};

#endif
