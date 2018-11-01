#ifndef DEFAULTPETSCFEATURE_SET
#define DEFAULTPETSCFEATURE_SET

#ifdef WITH_PETSCUI 

#include "typedefs.h"
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "petscdmda.h"

namespace SolverSelecter
{


class PetscTestingSpace
{
public:
  
    PetscTestingSpace( );

    virtual void
    extract_features(KSP &ksp,
                     std::map<std::string, double> &fmap,
                     int edge,
                     int internal,
                     bool matvec ) = 0;

    virtual void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) = 0;

    virtual void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double > &mmap) = 0;

    virtual void
    classify_measurements(std::map<std::string, double> &cmap ) = 0;

};


class DefaultPetscTestingSpace : public PetscTestingSpace
{
public:

    std::chrono::time_point<std::chrono::high_resolution_clock> time_start;
    PetscLogDouble space;


    void
    extract_features( KSP &ksp,
                      std::map<std::string, double> &fmap,
                      int edge,
                      int interior,
                      bool matvecs ) override;
    void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) override;

    void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double> &mmap ) override;

    void
    classify_measurements(std::map<std::string, double> &cmap ) override;



};

}
#endif

#endif
