
#ifndef _ML_WAFFLES_H
#define _ML_WAFFLES_H

/** \file MLWaff.h
 * \brief Implimentations for the ML waffles class. */
#include "MachineLearning.h"


/**
 * This class is a selector class for using a binary machine learning model in solver selection. This
 * will be used to load up a model that was pre-built externally. The binary format still needs to be
 * designed, but this shell should allow straight integration once that is done. 
 **/
class _ML_Binary : public _SS_MachineLearning
{
public:

    std::string filename;
    bool trained;

    /** Default constructer */
    _ML_Binary( std::string _filename) : _SS_MachineLearning(), filename(_filename), trained(false)
    {
    }


    /** Default destructor */
    ~_ML_Binary() {}

    _SS_ErrorFlag ResetModel( std::string _filename ) {
         
         if ( filename != _filename ) {
              filename = _filename;
              trained = false;
         }
         return _SS_error_flag;    
    }

    /**
     * Training is done externally in this system, so instead, just load the model. 
     **/
    _SS_ErrorFlag TrainSystem(  ) override { 
        
        loadModel();
        return _SS_error_flag;
    }
    
    /**
     * Classify a matrix based on its _SS_Features and return chosen solver
     **/
    _SS_ErrorFlag Classify( _SS_features_map &features /**< the feature set of the matrix */,
                            _SS_Solver &solver /**< output, a "good" solver for the problem */) override {
       loadModel() 
       return _SS_error_flag;
    }


    _SS_ErrorFlag loadModel() { 
      
      if ( ! trained ) {
          
          trained = true;      
      }
      return _SS_error_flag;
    }

};

#endif

