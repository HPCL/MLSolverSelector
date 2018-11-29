

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
    raplMeasurement.reset();
    raplMeasurement.startPowerMonitoring();
}

void
DefaultPetscTestingSpace::stop_measurements(KSP &ksp,
        Vec &x,
        Vec &b,
        std::map<std::string, double> &mmap )
{
    auto now = std::chrono::high_resolution_clock::now();
    auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(now-time_start).count();
    
    raplMeasurement.stopPowerMonitoring();

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    
    mmap["CPUTime"] =  ( reason > 0 ) ? (double) durr : -1.0 ;
    mmap["Power"] = 11; //raplMeasurement.getMaximumPowerConsumption();

}

void
DefaultPetscTestingSpace::classify_measurements( std::map<std::string, double> &cmap )
{
    cmap["CPUTime"] = 0.3;
}

}
#endif 
