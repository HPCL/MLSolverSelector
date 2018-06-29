#ifndef SS_MACHINELEARNING_H
#define SS_MACHINELEARNING_H

#include "DatabaseInterface.h"

namespace SolverSelecter { 

  class MachineLearningInterface : public SSBase 
  {
  public:
      
      std::string TAG = "MachineLearning";
      std::map< std::string , std::pair<std::vector<std::string>, std::string> > parameters;
      std::shared_ptr<DatabaseInterface> database; /**< pointer to the database */

  public:

      /**
       * Constructor
       **/
      MachineLearningInterface();
      
      virtual ~MachineLearningInterface();

      virtual ErrorFlag BuildModel() = 0;      

      virtual ErrorFlag Serialize(std::string output) = 0;

      virtual ErrorFlag ClassifyImpl( features_map &features, Solver &solver ) = 0;
      
      virtual ErrorFlag CrossValidate(std::string algorithm, 
                                      std::vector< std::string > &features ) = 0;
      
      virtual ErrorFlag Classify( features_map &features, Solver &solver );
   

      virtual ErrorFlag 
      GetMachineLearningData(std::vector< int > &row_ids , 
                             std::vector< int > &solvers_labels, 
                             std::vector< std::string > &features_labels ,
                             std::vector< std::string > &classification_labels,
                             std::vector< int > &solvers_data,  
                             std::vector< std::vector< double > > &feature_data,
                            std::vector< std::vector< bool > > &classification_data );

      virtual ErrorFlag Initialize( std::shared_ptr< DatabaseInterface > _database);

      virtual ErrorFlag CrossValidateAll( std::vector< std::string > algorithms, bool all  ); 

      virtual ErrorFlag CVFeaturesSpace (std::string alg, int level , 
                                         std::vector<std::string> &fnames,  
                                         std::vector<std::string> &sub_features);
   
  };
}
#endif
