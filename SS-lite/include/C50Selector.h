#ifndef C50_SELCETOR_H
#define C50_SELCETOR_H  

#include <map>
#include <string>

#include "C50Predictor.h"
#include "Solvers.h"
#include <memory>

class C50Interface 
{
  public:
      
      std::unique_ptr<C50Predictor> predictor;
      std::map< int, std::string > solverMap; 
      std::string m_filestem;
      

      C50Interface( std::string filestem );       
      virtual ~C50Interface();

      // Load the model as defined in the file stem. 
      int loadModel();


      int StringSplit(const std::string &s, const char *delim, std::vector< std::string > &result );

      // Based on this feature set, return a good solver .  
      Solver Predict( std::map<std::string, double> &features );  
};

#endif
