#ifndef SS_MACHINELEARNING_H
#define SS_MACHINELEARNING_H

#include "ParameterBase.h"
#include "DataBase.h"

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
class _SS_MachineLearning : public _SS_Parameters
{
protected:

    bool trained;  /**< has the system been trained already? */
    std::shared_ptr<_SS_DataBaseBase> database; /**< pointer to the database */

public:

    /**
     * Constructor
     **/
    _SS_MachineLearning(  );

    /**
     * Destructor
     **/
    virtual ~_SS_MachineLearning();

    _SS_ErrorFlag Initialize( std::shared_ptr< _SS_DataBaseBase > _database );

    /** This must be implemented by the subclass to actually train the system */
    virtual _SS_ErrorFlag TrainSystem( std::vector< std::pair< std::string, std::string > > &sset /**< set of allowed solvers*/) = 0;

    /** This function should predict if the solver will be good or bad */
    virtual _SS_ErrorFlag Predict( _SS_Features &features, /**< feature set for the matrix */
                                   const int &hash /**< hash of the solver to test for */,
                                   std::vector<bool> &good /**< vector saying if the solver is good FIXME why vector? */ ) = 0;


    /**
     * Train the system based on a list. If sset is empty, all solvers in the database are fair game. If
     * the list is populated, only the solvers listed in sset can be used.
     **/
    _SS_ErrorFlag Train( std::vector < std::pair < std::string, std::string >> &sset /** set of solvers to use in training (empty for all) */);

    /**
     * Classify a matrix based on its _SS_Features and return chosen solver
     **/
    _SS_ErrorFlag Classify( _SS_Features &features /**< the feature set of the matrix */,
                            _SS_Solver &solver /**< output, a (hopefully) "good" solver for the problem */);

};

#endif
