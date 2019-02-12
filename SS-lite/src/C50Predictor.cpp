
#include "C50Predictor.h"
#include <iostream>

C50Predictor::C50Predictor(String filestem ) {
    cc.RULES = false ;
    cc.RULESUSED = false;
    cc.FileStem = filestem;
    cc.loadModel();

    for ( int i = 0; i < cc.MaxAtt+1; i++ ) {
        if (cc.AttName[i] && strcmp(cc.AttName[i], "N/A") ) {
            attributes.push_back(cc.AttName[i]);
            if ( attributes.back() == "solver" ) {
                for ( int j = 0; j < cc.MaxAttVal[i]; j++ ) {
                    if ( cc.AttValName[i][j] && strcmp(cc.AttValName[i][j] , "N/A")  ) {
                        solvers.push_back( std::stoi(cc.AttValName[i][j]) );
                    }
                }

            }
            if ( attributes.back() == "class" ) {
                classAttr = attributes.size()-1;
            }
        }

    }
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
    DataRec predict1 = cc.GetDataRecFromVec(vfeatures, false);
    return Predict(predict1);
    
}


int C50Predictor::findAGoodSolver(std::map<std::string, double> &features) {
   
      bool good = false; 

      std::map<std::string, std::string> sfeatures;
      for ( auto &it : features ) {
        sfeatures.insert(std::make_pair( it.first, std::to_string(it.second) ));
      }

      for ( auto &it : solvers ) {
          sfeatures["solver"] = std::to_string(it);         
          if ( Predict(sfeatures) ) {

            std::cout << "\t Tried Solver " << it << " and... \\(ʘ‿ʘ)/  \n " ; 
            return it ;
          }
          else 
            std::cout << "\t Tried Solver " << it << " and... ¯\\_(ツ)_/¯  \n " ; 
      } 
      return -1;      
}

bool C50Predictor::Predict(DataRec Case) {
    return std::strcmp( cc.ClassName[cc.Classify(Case,cc.GCEnv)] , "bad" );  
}