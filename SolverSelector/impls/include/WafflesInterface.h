#ifndef SS_WAFFLES_H
#define SS_WAFFLES_H

#include "MachineLearningInterface.h"

#if WITH_WAFFLES

#include "GClasses/GMatrix.h"
#include "GClasses/GLearner.h"
#include "GClasses/GKNN.h"
#include "GClasses/GDecisionTree.h"
#include "GClasses/GNaiveBayes.h"
#include "GClasses/GTransform.h"
#include "GClasses/GDom.h"
namespace SolverSelecter
{

class WafflesInterface : public MachineLearningInterface
{
private:

    GClasses::GSupervisedLearner *model;     /**< the machine learning model */
    GClasses::GMatrix            *features;  /**< waffles matrix containing the features */
    GClasses::GMatrix            *labels ;   /**< waffles matrix containing the labels */

    std::vector< int > solver_hash_list; 
    std::map< int, std::string > solver_hash_map;
    std::vector< std::string > features_order ; 
    std::vector< std::string > labels_order; 

    /** Default constructer */
public:
    WafflesInterface();
    ~WafflesInterface();
private:
    /** Serialze the model */
    ErrorFlag SerializeImpl(std::vector <std::string>&CNames,std::string output) override;
    ErrorFlag BuildImpl(std::vector <std::string>&CNames) override;
    ErrorFlag BuildImplFromFile(std::string &Filename) override;
    ErrorFlag ClassifyImpl( features_map &afeatures, Solver &solver) override;
    ErrorFlag CrossValidateImpl( int folds ) override;
   
    /* Build the model from the serialized version */ 
    ErrorFlag BuildModelFromSerial(std::string input);

    /** Build the model from the database at runtime */
    ErrorFlag BuildModelAtRuntime(std::vector <std::string>&CNames );

    /** Extra function that converts the data into a GMatrix */
    ErrorFlag ImportData(std::vector <std::string>&CNames, std::vector < std::string > &feature_list,
                         std::shared_ptr<GClasses::GMatrix> &matrix,
                         int &num_labels );

    /** Callback function for cross validation tests. */
    static void CrossValidateCallback(void* pSupLearner, size_t nRep, size_t nFold, double foldSSE, size_t rows);

};

}
#endif
#endif
