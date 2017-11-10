
#include "MachineLearning.h"

/**
 * \file MachineLearning.C
 * \brief Impl files for machine learning classes
 **/


_SS_MachineLearning::_SS_MachineLearning() : trained(false)
{
}

_SS_MachineLearning::~_SS_MachineLearning()
{
}

_SS_ErrorFlag _SS_MachineLearning::Initialize( std::shared_ptr< _SS_DataBaseBase > _database )
{
    database = _database;
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_MachineLearning::Train( std::vector < std::pair < std::string, std::string >> &sset /** set of solvers to use in training (empty for all) */)
{

    if ( !trained )
    {
        TrainSystem( sset );
    }
    trained = true;
    return _SS_error_flag;
}


_SS_ErrorFlag _SS_MachineLearning::Classify( _SS_Features &features /**< the feature set of the matrix */,
                                             _SS_Solver &solver /**< output, a (hopefully) "good" solver for the problem */)

{

    /* Classify the system, choose the best solver, get the parameter string for
     * the solver, and initialize the solver with it */
    std::vector< bool > good;

    bool found_one = false;
    std::vector< int > hash_list;
    database->GetSolverHashList(hash_list);

    for ( auto hash : hash_list )
    {
        Predict( features, hash, good );
        auto it = good.begin();
        while ( it != good.end() )
        {
            if (*it) it++;
            else break;
        }
        if ( it == good.end() )
        {
            found_one = true;
            database->GetSolver( hash, solver );

            break;
        }
    }
    if ( !found_one ) std::cout<<"No Good Solvers -- YOYO\n";

    return _SS_error_flag;

}

