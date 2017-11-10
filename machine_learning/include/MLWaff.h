#ifndef _ML_WAFFLES_H
#define _ML_WAFFLES_H

#ifdef WITH_WAFFLES

/** \file MLWaff.h 
 * \brief Implimentations for the ML waffles class. */ 
#include "MachineLearning.h"
#include "ParameterBase.h"
#include "DataBase.h"
#include "Features.h"
#include "GClasses/GMatrix.h"
#include "GClasses/GLearner.h"
#include "GClasses/GKNN.h"
#include "GClasses/GDecisionTree.h"
#include "GClasses/GNaiveBayes.h"
#include "GClasses/GTransform.h"

namespace MachineLearning 
{

/**
 * The _ML_Waffles class does the actual ML for the Machine learning class. Need
 * to change the name to _SS_Waffles for consistency 
 **/
class _ML_Waffles : public _SS_MachineLearning  
{
  public:

      std::string method; /**< name of the method */
      GClasses::GSupervisedLearner *model;    /**< the machine learning model */
      GClasses::GAutoFilter        *fmodel;   /**< a filtered machine learning model */    
      GClasses::GMatrix            *features;  /**< waffles matrix containing the features */
      GClasses::GMatrix            *labels ;  /**< waffles matrix containing the labels */
      
      std::vector< std::string > attrs; /**< list of the attributes associated with the features matrix */
      
      int features_c; /**< number of feature columns */
      int features_r; /**< number of feature rows */
      int labels_r;   /**< number of label rows */
      int labels_c;   /**< number of label columns */

      /** Default constructer */
      _ML_Waffles();

      /** Default destructor */
     ~_ML_Waffles();

     /**
      * Train the chosen model. The training set can be restricted to only certain solver/preconditioner 
      * pairs. If sset.size() > 0, only this pairs listed in sset are used for training, and hence, for 
      * solving the problem. Future versions might allow for restriction based on parameters as well.
      **/
      _SS_ErrorFlag TrainSystem( std::vector< std::pair< std::string,std::string>> &sset 
                                 /**< optional set of solver/precond pairs to use */) override ;
      
      /**
       * Predict if a feature set is going to be a good solver. 
       **/ 
      _SS_ErrorFlag Predict( _SS_Features &features, /**< feature set to test */
                             const int &solver_hash,    /**< hash of the solver to test */
                             std::vector<bool> &good /**< bool staing if solver is good in each category */
                           ) override ;  

      /** Extra function that converts the data into a GMatrix */ 
      _SS_ErrorFlag ImportData( std::shared_ptr<GClasses::GMatrix> &matrix,
                         std::vector< int > &labelsdim,  
                         std::vector< std::pair< std::string, std::string > > &sset );
 

};
}


#endif
#endif

