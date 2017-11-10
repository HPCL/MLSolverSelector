
#ifndef SS_DATABASEBASE_H
#define SS_DATABASEBASE_H

/* API specific includes */
#include "Solvers.h"
#include "Measurements.h"
#include "Features.h"

/* Library files */
/**
 * \file DataBase.h
 * \brief Header files for integration with SQLITE3 
 **/

/**
 *  Class for reading and writing the database with SQLITE3. At the moment, only sqlite3 can be used
 *  for the database. In the future, might add support for other database types, so, this is written 
 *  as an abstact base class with a sqlite3 subclass. The base class contains pure virtual functions that
 *  are called at any point from other points in the code. **/
class _SS_DataBaseBase
{

public:
  
  std::string database_name; /**< filename of the database */

  /**
   *  Constructor. 
   **/
  _SS_DataBaseBase(const std::string &name /**< filename of the database */);

  /**
  *  Destructor 
  **/
  virtual ~_SS_DataBaseBase( );
  
   /**
   * Initialize the database. Initialization is completed outside the constructor to 
   * allow for some runtime options to be set. 
   **/
  virtual _SS_ErrorFlag Initialize() = 0 ; 
 
  /**
   * Finalize the database. 
   **/
  virtual _SS_ErrorFlag Finalize() = 0 ; 
    
  /**
   * WriteToFile is called at the end of the database building process. This is 
   *    
   **/
  virtual _SS_ErrorFlag WriteToFile() = 0;
  
  /**
   *  Add a row to the database based on a solver, matric and set of features and measurements 
   **/
  virtual _SS_ErrorFlag AddRow( const _SS_Solver &solver  /**< the solver used (includes parameters) */, 
                                const std::string &matrix  /**< the name of the matrix solved */,
                                const _SS_Features &features /**< the feature set extracted from the matrix */,
                                const _SS_Measurements &measurement /**< the measurements made during the solve */ ) = 0;
  
  /**
   * Classify the solvers in the database as good or bad. The measurements structure contains all the 
   * bad values required to do this classification.  
   **/
  virtual _SS_ErrorFlag ClassifySolvers( const _SS_Measurements &measure /**< input, used to get "bad" values */ ) = 0;
  
  /**
   * Retrive the solver from the database with the hash "shash" 
   **/
  virtual _SS_ErrorFlag GetSolver( int &shash /**< hash of the solver we want to get */, 
                                   _SS_Solver &solver /**< output, the solver set */ ) = 0;
  
   /**
   *  Returns a list of "shash" (i.e, solvers) in the database. This is used during machine learning
   *  to find a solver that is good.  
   **/
  virtual _SS_ErrorFlag GetSolverHashList( std::vector< int > &shashlist /**< output, vector of solvers in the database*/) = 0;

  /**
   *  Import the data from the database into matrix format for solving with the ML package. 
   **/
  virtual _SS_ErrorFlag ImportData( std::vector< std::pair< std::string , std::string > > & sset, 
                             std::vector< std::vector < std::string > > &data , 
                             std::vector<std::string> &column_names,
                             std::vector<int> &feature_or_label ) =0;

protected:
  
  /**
   *  These functions are implimented by the base class. These should not be changed in subclasses. This ensures
   *  that solver/preconditioners will have consistent names across all database implimentations. Maybe this isn't
   *  the best way to go about this. StableHashPstring takes a string and creates a hash. It is "Stable" because 
   *  the hash for a string is constant. Collisions shouldn't be a problem for the database sizes we are looking at
   **/
  _SS_ErrorFlag StableHashPString( const std::string &pstring /**< string to hash */, 
                                   int &phash /**< output, the hash, a positive integer */ );

  /**
   * Return the hash for the solver and matrix combination. Two hashes are returned. The uhash represents a 
   * hash for the solver applied to the matrix. This should be used as the rowid for the database. The shash 
   * is a hash for a particular solver. This hash can be used to search the database for all rows that use
   * a particular solver. 
   **/
  _SS_ErrorFlag GetHash( const _SS_Solver &solver, const std::string &matrix, int &uhash, int &shash );
};


#endif 



