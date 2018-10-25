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
MachineLearningInterface::Serialize(std::string outputfile) {
    SerializeImpl(outputfile);
    return 0;
}

ErrorFlag 
MachineLearningInterface::Build() {

    if ( !modelBuilt) {
        BuildImpl();
        modelBuilt = true;
    }
    return 0;
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
   Build();
   if ( features.size() > 0 )
      ClassifyImpl(features,solver);
   else 
      solver.Clear();
   return 0;
}

ErrorFlag MachineLearningInterface::AddToBanedListAndTryAgain( features_map &afeatures,
                                                               Solver &solver ) 
{
  Build();
  int solverHash;
  HashString(solver.solver, solverHash);
  bannedList.insert(solverHash);
  solver.Clear();
  Classify(afeatures, solver);
  return 0;
}

bool MachineLearningInterface::isBaned(Solver &solver) {
  int solverHash;
  HashString(solver.solver, solverHash);
  return ( bannedList.find(solverHash) != bannedList.end() );
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
  
  if ( features_labels.size() == 0 || row_ids.size() == 0 
                                   || classification_labels.size() == 0
                                   || solvers_labels.size() == 0 )
  {
      throw SSException("Database Failure"); 
  }

  return error_flag;
}

ErrorFlag
MachineLearningInterface::CrossValidate( int folds  )
{
   Build(); 
   CrossValidateImpl( folds );
   return 0;
}


}
