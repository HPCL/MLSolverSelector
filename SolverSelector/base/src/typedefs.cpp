#include "typedefs.h"

namespace SolverSelecter
{

ErrorFlag
StringSplit(const std::string &s,
            const char *delim,
            std::vector< std::string > &result )
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while ( std::getline( ss,item,delim[0]) )
    {
        if (!item.empty())
            result.push_back( item );
    }
    return error_flag;
}

ErrorFlag
HashString(const std::string &pstring,
           int &phash)
{
    const char * str = pstring.data();
    unsigned int h;
    unsigned char *p;
    unsigned int MULTIPLIER = 31;
    h = 0;
    for (p = (unsigned char*)str; *p != '\0'; p++)
        h = MULTIPLIER * h + *p;
    phash = abs( (int) h ) ;
    return error_flag;
}
}

