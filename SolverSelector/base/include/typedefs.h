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

namespace SolverSelecter { 

/* Put these here now, so they are easy to change/overide if needed */
typedef unsigned int ErrorFlag ;
static ErrorFlag error_success = 0;
static ErrorFlag error_fail = 1;
static ErrorFlag error_flag = 0; // TODO 

typedef std::map< std::string, std::set< std::string > > parameters_map; 
typedef std::map< std::string, double > features_map;
typedef std::map< std::string, double > measurements_map;

struct SSException : public std::exception 
{
  std::string s;
  SSException(std::string ss) : s(ss) {} 
  ~SSException() {}
  const char* what() const throw() { return s.c_str(); } 
};

class SSBase 
{

  public:  
  std::string TAG; /// Tag to associate class with a ascii input variable
  std::map< std::string, std::string > parameter_help, parameter_values;

  SSBase(std::string tag) : TAG(tag) {} 

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
    if ( valid == 1 )
        parameter_values[name] = value;
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

ErrorFlag 
StringSplit(const std::string &s,
            const char *delim,
            std::vector< std::string > &result );

ErrorFlag
HashString( const std::string &pstring,
                              int &phash );


}

#endif

