#ifndef PETSCINTERFACE_H
#define PETSCINTERFACE_H

#include <string>
#include <map>
#include <memory>

#include "Solvers.h"
#include "C50Selector.h"
#include "petscksp.h"
#include "petsc/private/kspimpl.h"

// This class provides a selection method for chosing a different feature set 

class PetscInterface 
{
public:
    
    KSP _ksp; 
    std::string internal_prefix;
  
    Solver currentSolver; 

    std::string filestem; // filestem for the C50 model.  
    std::string featureSet; // Name of the feature set to use. 

    int edgeSamples = 20;
    int interiorSamples = 40;
    bool useMatVecs = true;  // Should we use matvecs (true) or memory location (false) 
      
    std::unique_ptr<C50Interface> c50Interface;

    // to extract the sample columns from the matrix. 
    PetscInterface();

    virtual int Initialize(std::string filestem_, std::string featureSet_, int edge, int interior, bool matvec);

    virtual int ClassifyAndSolve(KSP &ksp);

    virtual int SolveSystem( KSP &ksp, Solver &solver, bool &success);
    
    virtual int SetSolver( KSP ksp, Solver solver ); 

    virtual int ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap) ;
        
    virtual int GetDefaultSolver( Solver &solver ) ;
   
    virtual int Predict(std::map<std::string, double> &features_map, Solver &solver);

    virtual int PrintSolver();

    virtual  ~PetscInterface();
  
private: 
    void PetscCopyFunction( KSP nksp, KSP oksp );   

};

class PetscCoupler {
  public : 
      static bool initialized;  
      static std::shared_ptr<PetscInterface> static_petcs_interface_ptr; 
    
      static PetscErrorCode monitor_ss_convergence(KSP ksp, PetscInt iter, PetscReal rnorm , void * ctx ) ;
      static PetscErrorCode KSPSetUp_SS(KSP ksp);
      static PetscErrorCode KSPSolve_SS(KSP ksp);
      static PetscErrorCode KSPDestroy_SS(KSP ksp);
      static PetscErrorCode KSPView_SS(KSP ksp, PetscViewer viewer);
      static PetscErrorCode KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp);
      static PetscErrorCode KSPCreate_SS(KSP ksp);

      static void CreateSolverSelectorKSP() {
          KSPRegister("KSPSS", KSPCreate_SS);
      }

};

#endif

