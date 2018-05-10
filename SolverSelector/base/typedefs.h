#ifndef SS_TYPEDEF_H
#define SS_TYPEDEF_H

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

/* Put these here now, so they are easy to change/overide if needed */
typedef unsigned int _SS_ErrorFlag ;
static _SS_ErrorFlag _SS_error_flag = 0;
static _SS_ErrorFlag _SS_error_fail = 1;

typedef std::map< std::string, std::set< std::string > > _SS_parameters_map; 
typedef std::map< std::string, double > _SS_features_map;
typedef std::map< std::string, double > _SS_measurements_map;


namespace _SS_Utils { 

/**
 * Split the strings based on a delimiter
**/
_SS_ErrorFlag StringSplit( const std::string &s, const char *delim, std::vector< std::string > &result ) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while ( std::getline( ss,item,delim[0]) )
        {
            if (!item.empty())
                result.push_back( item );
        }
        return _SS_error_flag;
};

    /**
     * Get a stable hash for a string  
     ***/
    _SS_ErrorFlag HashString( const std::string &pstring /**< string to hash */,
                                     int &phash /**< output, the hash, a positive integer */ )
    {
        const char * str = pstring.data();
        unsigned int h;
        unsigned char *p;
        unsigned int MULTIPLIER = 31;
        h = 0;
        for (p = (unsigned char*)str; *p != '\0'; p++)
            h = MULTIPLIER * h + *p;
        phash = abs( (int) h ) ;
        return _SS_error_flag;
    }


}


#endif

