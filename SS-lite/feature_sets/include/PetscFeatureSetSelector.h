
#ifndef PETSCFEATURESETSELECTOR_H
#define PETSCFEATURESETSELECTOR_H

//# Feature Set Include 
#include "example_featureExtraction.h" 

// Feature Set Include 
#include "full_featureExtraction.h" 

//__SEARCH__INCLUDE

// This is file is automatically modified whenever a new feature set it added using the generate
// command. 

namespace PetscFeatureSetSelector {

int ExtractJacobianFeatures( std::string featureSet, Mat J , int edge, int interior, std::map<  std::string, double > &fnames , bool matvecs=true )
{

    //Feature Set 
    if ( featureSet == "example" ) return FeatureExtraction_example::ExtractJacobianFeatures(J,edge,interior,fnames,matvecs); 

    //Feature Set 
    if ( featureSet == "full" ) return FeatureExtraction_full::ExtractJacobianFeatures(J,edge,interior,fnames,matvecs); 

    //__SEARCH__IF ; 
    
    std::cout << " Could not find a feature set with the name " << featureSet << ". Aborting " ;
    std::abort();
  
}

bool featureSetExists(std::string featureSet) {
   
//   Feature Set Exists 
   if( featureSet == "example") return true; 

//   Feature Set Exists 
   if( featureSet == "full") return true; 

//__SEARCH__EXISTS 
    return false;
}



}

#endif

