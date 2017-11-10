
#include "DataBase.h"

/**
 * \file DataBase.C
 * \brief Implimentation files for database support
 **/

_SS_DataBaseBase::_SS_DataBaseBase(const std::string &name /**< filename of the database */)
{
    database_name = name;
}

_SS_DataBaseBase::~_SS_DataBaseBase( )
{

}


_SS_ErrorFlag _SS_DataBaseBase::StableHashPString( const std::string &pstring /**< string to hash */,
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

_SS_ErrorFlag _SS_DataBaseBase::GetHash( const _SS_Solver &solver, 
                                         const std::string &matrix, 
                                         int &uhash, 
                                         int &shash )
{
    std::string solverstring;
    solver.GetSolverString( solverstring );
    StableHashPString( solverstring, shash );
    solverstring = matrix + " " + solverstring;
    StableHashPString( solverstring, uhash );
    return _SS_error_flag;
}



