#include "MachineLearningInterface.h"

namespace SolverSelecter
{

MachineLearningInterface::MachineLearningInterface() : SSBase("MACHINELEARNING")
{
}

MachineLearningInterface::~MachineLearningInterface()
{

}


ErrorFlag
MachineLearningInterface::Initialize( std::shared_ptr< DatabaseInterface > _database)
{
    database = _database;
    return error_flag;
}
    
ErrorFlag 
MachineLearningInterface::Classify( features_map &features,
                                    Solver &solver ) 
{
   if ( features.size() > 0 )
      ClassifyImpl(features,solver);
   else 
      solver.Clear();
   return 0;
}

ErrorFlag
MachineLearningInterface::GetMachineLearningData(std::vector< int > &row_ids , 
                                                 std::vector< int > &solvers_labels, 
                                                 std::vector< std::string > &features_labels ,
                                                 std::vector< std::string > &classification_labels,
                                                 std::vector< int > &solvers_data,  
                                                 std::vector< std::vector< double > > &feature_data,
                                                 std::vector< std::vector< bool > > &classification_data ) 
{
  database->GetMachineLearningData(row_ids,solvers_labels,features_labels, classification_labels,
                                   solvers_data, feature_data, classification_data );
  
  if ( features_labels.size() == 0 || row_ids.size() == 0 || classification_labels.size() == 0 || solvers_labels.size() == 0 )
  {
      std::cout << " Database Import Failure: Default solver used from this point on \n";
      throw SSException("Database Failure"); 
  }

  return error_flag;
}

ErrorFlag
MachineLearningInterface::CrossValidateAll( std::vector< std::string > algorithms, bool all  )
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

    return error_flag;

}

ErrorFlag
MachineLearningInterface::CVFeaturesSpace ( std::string alg,
        int level ,
        std::vector<std::string> &fnames,
        std::vector<std::string> &sub_features)
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
    return error_flag;
}

}
