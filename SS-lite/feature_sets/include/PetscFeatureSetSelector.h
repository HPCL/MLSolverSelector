
#ifndef PETSCFEATURESETSELECTOR_H
#define PETSCFEATURESETSELECTOR_H
#include "petscdmda.h"
#include "petscksp.h"
#include <iostream>


// Feature Set Include 
 namespace FeatureExtraction_example { int ExtractJacobianFeatures(Mat J, int edge, int interior, std::map<std::string, double> &fnames, bool usePercent, bool matvecs); }  

// Feature Set Include 
 namespace FeatureExtraction_full { int ExtractJacobianFeatures(Mat J, int edge, int interior, std::map<std::string, double> &fnames, bool usePercent, bool matvecs); }  

//__SEARCH__INCLUDE


namespace PetscFeatureSetSelector {

  int ExtractJacobianFeatures( std::string featureSet, Mat J , int edge, int interior, std::map<  std::string, double > &fnames , bool usePercent = false, bool matvecs=true )
  { 

      //Feature Set 
    if ( featureSet == "example" ) return FeatureExtraction_example::ExtractJacobianFeatures(J,edge,interior,fnames,usePercent, matvecs); 

    //Feature Set 
    if ( featureSet == "full" ) return FeatureExtraction_full::ExtractJacobianFeatures(J,edge,interior,fnames,usePercent, matvecs); 

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

