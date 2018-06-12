#ifndef SS_MACHINELEARNING_H
#define SS_MACHINELEARNING_H

#include "DatabaseInterface.h"

/**
 * \file MachineLearning.h
 * \brief Header files for machine learning classes
 **/


/******************************************************************************
 **************************MACHINE LEARNING CLASS*****************************/

/**
 * This class acts as a middle man between the solver and the machine learning code. The _ML_Waff
 * class is currently used to complete that actual machine learning. This is an abstract class
 * that must be implimented by the actual ML code.  **/
class _SS_MachineLearning
{
public:

    bool trained;  /**< has the system been trained already? */
    std::shared_ptr<_SS_DatabaseInterface> database; /**< pointer to the database */

public:

    /**
     * Constructor
     **/
    _SS_MachineLearning(  ) : trained(false)
    {
    }

    /**
     * Destructor
     **/
    virtual ~_SS_MachineLearning()
    {
    }

    _SS_ErrorFlag Initialize( std::shared_ptr< _SS_DatabaseInterface > _database )
    {
        database = _database;
        return _SS_error_flag;
    }

    /** This must be implemented by the subclass to actually train the system */
    virtual _SS_ErrorFlag TrainSystem( ) = 0;

    /**
     * Classify a matrix based on its _SS_Features and return chosen solver
     **/
    virtual _SS_ErrorFlag Classify( _SS_features_map &features /**< the feature set of the matrix */,
                            _SS_Solver &solver /**< output, a (hopefully) "good" solver for the problem */) = 0;

    
    virtual _SS_ErrorFlag CrossValidate( std::string algorithm, std::vector< std::string > &features )
    {
        return _SS_error_flag;
    }

    virtual _SS_ErrorFlag CrossValidateAll( std::vector< std::string > algorithms, bool all  )
    {
        for ( auto it : algorithms )
        {
            if ( all )
            {
                /* Recurse over all posible feature combos */
                std::vector< std::string > fnames, subfeatures ;
                CVFeaturesSpace( it, 0, fnames, subfeatures );
            }
            else
            {
                /* just use the main feature set */
                std::vector< std::string > fnames ;
                database->GetFeatureLabels(fnames);
                CrossValidate( it, fnames );
            }
        }
        return _SS_error_flag;

    }

    virtual _SS_ErrorFlag CVFeaturesSpace ( std::string alg, int level , std::vector<std::string> &fnames,  std::vector<std::string> &sub_features)
    {
        if (level == 0)
            database->GetFeatureLabels(fnames);

        while ( fnames.size() > 0 )
        {
            std::vector< std::string > fnames_copy = fnames;
            sub_features.push_back( fnames.back() ) ;
            fnames_copy.pop_back();

            CrossValidate( alg, sub_features );

            if ( fnames_copy.size() > 0 )
                CVFeaturesSpace( alg, level+1, fnames_copy, sub_features ) ;
            sub_features.pop_back();
            fnames.pop_back();
        }
        return _SS_error_flag;
    }
    /**
     * Train the system based on a list. If sset is empty, all solvers in the database are fair game. If
     * the list is populated, only the solvers listed in sset can be used.
     **/
    _SS_ErrorFlag Train( )
    {

        if ( !trained )
        {
            TrainSystem();
        }
        trained = true;
        return _SS_error_flag;
    }


    };

#endif
