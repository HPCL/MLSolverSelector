

#if WITH_PETSCUI 

#include "PetscFeatures.h"
#include "FeatureExtract.h"


namespace SolverSelecter
{

PetscTestingSpace::PetscTestingSpace() {}


#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "petscdmda.h"

DefaultPetscTestingSpace::DefaultPetscTestingSpace(std::chrono::nanoseconds interval) : PetscTestingSpace(), raplMeasurement(interval) {}

void
DefaultPetscTestingSpace::extract_features(KSP &ksp,
        std::map<std::string, double> &fmap, int edge, int internal, bool matvec )
{
    Mat AA,PP;
    KSPGetOperators(ksp, &AA, &PP );
    std::cout << edge << internal << matvec << std::endl; 
    ExtractJacobianFeatures(AA,edge,internal,fmap,matvec);
}

void
DefaultPetscTestingSpace::start_measurements(KSP &ksp,
        Vec &x,
        Vec &b )

{
    time_start = std::chrono::high_resolution_clock::now();
   

   
   //memorymeasure.startMemMonitoring();
  papimeasure.startPowerMonitoring();
  //  perfmeasure.startMonitoring();
   cachemeasure.startMemMonitoring();
  // perfmeasure.startMonitoring();
}

void
DefaultPetscTestingSpace::stop_measurements(KSP &ksp,
        Vec &x,
        Vec &b,
        std::map<std::string, double> &mmap )
{
    auto now = std::chrono::high_resolution_clock::now();
    auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(now-time_start).count();
    
   papimeasure.stopPowerMonitoring(mmap);
   //memorymeasure.stopMemMonitoring(mmap);
   cachemeasure.stopMemMonitoring(mmap);
  // perfmeasure.stopMonitoring(mmap); 

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    
   mmap["CPUTime"] = (double) durr;
    printf("Time  %.12f \n", (mmap["CPUTime"]) * 1e-9);
    if (reason <= 0)
    {
        for (auto it : mmap)
        {
            mmap[it.first] = -1;
        }
    }
   mmap["PACKAGE_ENERGY"] =  ( reason > 0 ) ? (double) durr : -1 ;
   mmap["DRAM_ENERGY"] =  ( reason > 0 ) ? (double) durr : -1 ;
   //mmap["CPUTime"] = mmap["PACKAGE_ENERGY"];
   // mmap["Power"] =10;
   // printf("pwer test %.2f\n", raplMeasurement.getpower() );
}

void
DefaultPetscTestingSpace::classify_measurements( std::map<std::string, double> &cmap )
{
   //cmap["CPUTime"] = 0.1;
   // cmap["PACKAGE_ENERGY"] = 0.1;
   // cmap["DRAM_ENERGY"] = 0.3;
    //cmap["CPUTime"] = 0.3;
}

}
#endif 
