#ifndef HEADER_DUMMY_FEATURE_SET
#define HEADER_DUMMY_FEATURE_SET

#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "MLWaff.h"

class DummyTestingSpace : public PetscTestingSpace 
{
public:
  std::chrono::time_point<std::chrono::high_resolution_clock> time_start;
  
   void set_machine_learning_model( std::shared_ptr<_SS_MachineLearning> &machine ) {
    
      machine.reset( new _ML_Waffles() ) ; 
      // machine.reset( new _ML_Binary( filename ) ; // not implemented yet. 
   } 

   void extract_features( KSP &ksp, std::map<std::string, double> &fmap ) override
  {      
      /* This is a dummy feature set. Basically sets some random junk */
      /* This is where we would insert the Petsc Matrix-free extraction algorithm. 
       */
       
      fmap["happy"] = (double) rand();
      fmap["sad"] = (double) -1.0*rand();
      fmap["mad"] = (double) rand() / RAND_MAX ; 
  }
 
  void start_measurements( KSP &ksp, Vec &x, Vec &b ) override
  {      
      time_start = std::chrono::high_resolution_clock::now();
  }

  void stop_measurements( KSP &ksp, Vec &x, Vec &b, std::map<std::string, double> &mmap ) override
  {
        auto now = std::chrono::high_resolution_clock::now();
        auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(now-time_start).count();
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        mmap["CPUTime"] =  ( reason > 0 ) ? (double) durr : 1e300 ;
  }
  void classify_measurements( std::map<std::string, double> &cmap ) override
  {
       cmap["CPUTime"] = 0.3; 
  }
    
};



#endif

