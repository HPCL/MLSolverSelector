#ifndef SS_TYPEDEF_H
#define SS_TYPEDEF_H

#include <cstring>
#include <string>
#include <set>
#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>
#include <memory>
#include <limits>
#include <chrono>
#include <mpi.h> // MPI Support 

/** This is just a general file defining some important classes. This file 
 * is some way or another included in every compilation block when building 
 * the library */

namespace SolverSelecter { 

/* Put these here now, so they are easy to change/overide if needed. Error handling
 * is a big TODO. For the most part, most functions return a error_flag, although it 
 * is almost never set false, and is never actually tested.  */
typedef unsigned int ErrorFlag ;
static ErrorFlag error_success = 0;
static ErrorFlag error_fail = 1;
static ErrorFlag error_flag = 0; 

/** some typedefs used throughout to save typeing */
typedef std::map< std::string, std::set< std::string > > parameters_map; 
typedef std::map< std::string, double > features_map;
typedef std::map< std::string, double > measurements_map;

/** Another attempt at error handling. This time, a simple Exception we can throw when something
 * we don't like happens. This is used every now and then, but needs to be used more */
struct SSException : public std::exception 
{
  std::string s;
  SSException(std::string ss) : s(ss) {} 
  ~SSException() {}
  const char* what() const throw() { return s.c_str(); } 
};

/** This is the main base class that handles all the parameters. Classes what want to 
 * have parameter sets override this class. */ 
class SSBase 
{

  public:  
  std::map< std::string, std::string > parameter_help, parameter_values;

  SSBase(std::string tag) {} 

  std::map<std::string, std::string > GetParameters() {
    return parameter_values;
  }

  std::string GetParameterHelp(const std::string &parameter ) {
    if ( parameter_help.find(parameter) == parameter_help.end() ) 
      return "unknown parameter";
    return parameter_help[parameter];
  }

  ErrorFlag DumpParameterHelpMessage() {
    std::cout << " The parameters available for this interface are:\n" ;
    for ( auto &it : parameter_values ) {
      std::cout << "\t" << it.first << ": < " << it.second << "> : " << parameter_help[it.first] << std::endl;
    }
    return 0;
  }


  ErrorFlag Parse( std::vector< std::pair< std::string, std::string > > &p )
  {
    for ( auto it : p ) 
      SetParameter(it.first,it.second);
    return error_flag;
  }
  
  ErrorFlag AddParameter(std::string name, std::string value, std::string help )
  {
      parameter_help.insert(std::make_pair(name,help));
      parameter_values.insert(std::make_pair(name, value) );
      return error_flag;
  }

  ErrorFlag CheckValid(std::string name, std::string value , int &valid )
  {    
    valid = ( parameter_values.find(name) == parameter_values.end() ) ? 0 : 1; 
    return error_flag;
  }

  ErrorFlag SetParameter(std::string name, std::string value) 
  {
   
    int valid;
    CheckValid(name, value, valid) ;
    if ( valid == 1 ) {
        parameter_values[name] = value;
    } else {
    }
   return error_flag; 
  }

  std::string GetParameter(std::string name )
  {
     if ( parameter_values.find(name) != parameter_values.end() ) 
        return parameter_values[name]; 
     else
        throw SSException("Parameter Not Found"); 
  }

  bool GetParameterAsBool(std::string name) 
  {
    std::string str = GetParameter(name);
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
  }

  int GetParameterAsInt(std::string name)
  {
    return std::atoi(GetParameter(name).c_str());
  }

  double GetParameterAsDouble(std::string name)
  {
      return std::atof(GetParameter(name).c_str());
  }

};

/** Utility funtion to split a string on a delimiter */
ErrorFlag 
StringSplit(const std::string &s,
            const char *delim,
            std::vector< std::string > &result );

/** Utiliity function to repoducibly hash a string */ 
ErrorFlag
HashString( const std::string &pstring,
                              int &phash );


}

#endif

