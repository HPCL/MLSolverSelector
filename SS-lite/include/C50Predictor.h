#ifndef SS_C50PREDICTOR_H
#define SS_C50PREDICTOR_H

#include <map>
#include <string>
#include <vector>
#include <string>
#include <cstring>

#include "C50Core.h" // Defines String and C50Classifier 

class C50Predictor {
public:
    C50Classifier cc;
    std::vector<std::string> attributes;
    std::vector<int> solvers;
    std::vector<int>::iterator solversIter;
    
    int classAttr = -1;
    
    C50Predictor(String filestem);     
    std::vector<std::string> MapToVec( std::map<std::string, std::string> &features ); 
    
    bool Predict(DataRec Case);
    bool Predict(std::map<std::string, std::string> &features);
   
    bool trySolver(int solver, std::map<std::string, std::string> &sfeatures);
 
    int findAGoodSolver(std::map<std::string, double> &features ); 

};


#endif
