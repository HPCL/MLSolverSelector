#ifndef SS_MACHINELEARNING_H
#define SS_MACHINELEARNING_H

#include "DatabaseInterface.h"

namespace SolverSelecter { 
  
  /** Base class for interacting with external machine learning algorithms and packages */
  class MachineLearningInterface : public SSBase 
  {
  protected:
      
      std::shared_ptr<DatabaseInterface> database; /**< pointer to the database */
      std::set<int> bannedList;      
      bool modelBuilt = false;

  public:
      /** Constructor */
      MachineLearningInterface();

      /** Desctuctor */
      virtual ~MachineLearningInterface();
      
      // Build the model then serialize it to file. 
      ErrorFlag Serialize(std::string outputfile);

      // Initialize the machine learning model (but don't build it)
      ErrorFlag Initialize( std::shared_ptr< DatabaseInterface > _database);
      
      // Classify the feature set and return a "good" solver. 
      ErrorFlag Classify( features_map &features, Solver &solver );      
      
      // If the returned solver is bad, call this function to add it to the ignore list
      // and try again
      ErrorFlag AddToBanedListAndTryAgain( features_map &afeatures,
                                                   Solver &solver)  ;
      
      // Cross validate the given machine learning model. 
      ErrorFlag CrossValidate( int folds ); 
      
  
  private:

      // Build the machine learning model. 
      ErrorFlag Build(); 
      
      /* Build the machine learning model implementation*/
      virtual ErrorFlag BuildImpl() = 0;
      
      /* Classify the given feature set and return a good solver */
      virtual ErrorFlag ClassifyImpl( features_map &features, Solver &solver ) = 0;
      
      /* Build the model for the given algorithm */
      virtual ErrorFlag CrossValidateImpl( int folds ) = 0;
      
      /** Implementation for Serializing the model */
      virtual ErrorFlag SerializeImpl(std::string outputfile) = 0;
  
  protected:    
     
      /** Get the machine learning data from the database. */ 
      ErrorFlag 
      GetMachineLearningData(std::vector< int > &row_ids , 
                             std::vector< int > &solvers_labels, 
                             std::vector< std::string > &features_labels ,
                             std::vector< std::string > &classification_labels,
                             std::vector< int > &solvers_data,  
                             std::vector< std::vector< double > > &feature_data,
                            std::vector< std::vector< bool > > &classification_data );
      
      /** Is the given solver on the banned list */
      bool isBaned(Solver &solver);  

  };
}
#endif
