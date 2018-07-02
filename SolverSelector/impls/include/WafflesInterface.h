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
public:

    GClasses::GSupervisedLearner *model;     /**< the machine learning model */
    GClasses::GMatrix            *features;  /**< waffles matrix containing the features */
    GClasses::GMatrix            *labels ;   /**< waffles matrix containing the labels */

    std::vector< int > solver_hash_list; 
    std::map< int, std::string > solver_hash_map;
    std::vector< std::string > features_order ; 
    std::vector< std::string > labels_order; 


    /** Default constructer */
    WafflesInterface();


    /** Default destructor */
    ~WafflesInterface();


    ErrorFlag Serialize(std::string output) override;

    /**
     * Train the chosen model. The training set can be restricted to only certain solver/preconditioner
     * pairs. If sset.size() > 0, only this pairs listed in sset are used for training, and hence, for
     * solving the problem. Future versions might allow for restriction based on parameters as well.
     **/
    ErrorFlag BuildModel();
    ErrorFlag BuildModelFromSerial(std::string input);

    /**
     * Classify a matrix based on its Features and return chosen solver
     **/
    ErrorFlag ClassifyImpl( features_map &afeatures /**< the feature set of the matrix */,
                        Solver &solver /**< output, a (hopefully) "good" solver for the problem */) override;


    /** Extra function that converts the data into a GMatrix */
    ErrorFlag ImportData(
        std::vector < std::string > &feature_list,
        std::shared_ptr<GClasses::GMatrix> &matrix,
        int &num_labels );

    ErrorFlag CrossValidate( std::string algorithm, std::vector< std::string > &features_list ) override;
  
    ErrorFlag BuildModelAtRuntime( std::string algorithm, std::vector<std::string> &fnames );
    
    static void CrossValidateCallback(void* pSupLearner, size_t nRep, size_t nFold, double foldSSE, size_t rows);

};

}
#endif
#endif
