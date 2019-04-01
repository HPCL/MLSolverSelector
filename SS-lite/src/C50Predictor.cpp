
#include "C50Predictor.h"
#include <iostream>

// Only inlcude the core here, so we can encapsulate all the garbage, and not 
// have it pulled into the executable when compiling. 
#include "C50Core.h"

C50Predictor::C50Predictor(std::string filestem ) {
    
    cc = new C50Classifier();
  
    cc->RULES = false ;
    cc->RULESUSED = false;
    cc->FileStem = filestem.c_str();
    cc->loadModel();

    for ( int i = 0; i < cc->MaxAtt+1; i++ ) {
        if (cc->AttName[i] && strcmp(cc->AttName[i], "N/A") ) {
            attributes.push_back(cc->AttName[i]);
            if ( attributes.back() == "solver" ) {
                for ( int j = 0; j < cc->MaxAttVal[i]; j++ ) {
                    if ( cc->AttValName[i][j] && strcmp(cc->AttValName[i][j] , "N/A")  ) {
                        solvers.push_back( std::stoi(cc->AttValName[i][j]) );
                    }
                }

            }
            if ( attributes.back() == "class" ) {
                classAttr = attributes.size()-1;
            }
        }

    }
    solversIter = solvers.begin();


}

C50Predictor::~C50Predictor() {
  delete cc;
}

std::vector<std::string> C50Predictor::MapToVec( std::map<std::string, std::string> &features ) {
    std::vector< std::string > vfeatures;
    for ( auto &it : attributes ) {
        if ( features.find(it) != features.end() ) {
            vfeatures.push_back( features[it] );
        } else
            vfeatures.push_back("?");
    }
    return vfeatures;
}

bool C50Predictor::Predict(std::map<std::string, std::string> &features) {
    std::vector<std::string> vfeatures = MapToVec(features);
    DataRec predict1 = cc->GetDataRecFromVec(vfeatures, false);
    return std::strcmp( cc->ClassName[cc->Classify(predict1,cc->GCEnv)] , "bad" );  
    
}


int C50Predictor::findAGoodSolver(std::map<std::string, double> &features) {
   
      bool good = false; 

      std::map<std::string, std::string> sfeatures;
      for ( auto &it : features ) {
        sfeatures.insert(std::make_pair( it.first, std::to_string(it.second) ));
      }
      
      auto startIter = solversIter; 
      bool start = true;      
      
      while( start || solversIter != startIter ) {
          start = false;
          if (solversIter == solvers.end() ) 
             solversIter == solvers.begin(); 
          else if ( trySolver( *solversIter, sfeatures ) ) 
              return *solversIter;
          else {
              solversIter++;
          }      
      }  
      std::cout << " \t Could not find a single good solver \n " ;     
      return -1;      
}

bool C50Predictor::trySolver(int solver, std::map<std::string, std::string> &sfeatures) {
    sfeatures["solver"] = std::to_string(solver);         
    if ( Predict(sfeatures) ) {
      std::cout << "\t Tried Solver " << *solversIter << " and... \\(ʘ‿ʘ)/  \n " ; 
      return true;
    }
    else {
      std::cout << "\t Tried Solver " << *solversIter << " and... ¯\\_(ツ)_/¯  \n " ; 
      return false;
    }
}

