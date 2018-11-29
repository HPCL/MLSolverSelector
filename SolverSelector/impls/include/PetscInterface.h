#ifndef PETSCINTERFACE_H
#define PETSCINTERFACE_H

#ifdef WITH_PETSCUI

#include "SolverSelecter.h"
#include "UserInterface.h"
#include "WafflesInterface.h"
#include "Sqlite3Interface.h"
#include "C50Interface.h"

#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "PetscFeatures.h"

namespace SolverSelecter {

typedef SolverSelecter<KSP,Vec> PetscSS;

class PetscUI : public UserInterface<KSP,Vec>
{
public:

    KSP _ksp; 
    std::shared_ptr< PetscTestingSpace > testing_space;
    std::string internal_prefix;
    bool use_internal_ksp = true;
    
    PetscUI();

    virtual ErrorFlag Initialize() override; 

    virtual ErrorFlag SetTestingSpace();

    virtual ErrorFlag GetDataBaseImplimentation(std::shared_ptr< DatabaseInterface > &database); 
     
    virtual ErrorFlag GetMachineLearningModel(std::shared_ptr< MachineLearningInterface > &machine);
    
    virtual ErrorFlag ChangeSolver( KSP &A, Vec &x, Vec &b , bool &change_it ); 

    virtual ErrorFlag SolveSystem( KSP &ksp, Vec &x, Vec &b, Solver &solver, bool &success) override;
    
    virtual ErrorFlag SetSolver( KSP ksp, Solver solver ); 

    virtual ErrorFlag ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap) override;
    
    virtual ErrorFlag StartMeasurements(KSP &A,  Vec &x,  Vec &b ) override;
    
    virtual ErrorFlag StopMeasurements( KSP &A, Vec &x, Vec &b, std::map<std::string, double> &mmap ) override; 

    virtual ErrorFlag ClassificationMeasurements( std::map< std::string, double > &cmap ) override; 

    virtual ErrorFlag InitMatrix( std::string filename, std::unique_ptr<KSP> &A ) override;
    
    virtual ErrorFlag InitVector( const KSP &A, std::unique_ptr<Vec> &x ) override;
    
    virtual ErrorFlag CopyVector( const Vec &x, std::unique_ptr<Vec> &xx );
    
    virtual ErrorFlag FreeVector( std::unique_ptr<Vec> &x) override;
    
    virtual ErrorFlag FreeMatrix( std::unique_ptr<KSP> &A ) override;
    
    virtual ErrorFlag SetVector( std::unique_ptr<Vec> &x, std::string type_ ) override;
    
    virtual ErrorFlag GetDefaultSolver( Solver &solver ) override;
   
    virtual int GetNNZ( KSP &A) override; 

    virtual  ~PetscUI();
  
private: 
    void PetscCopyFunction( KSP nksp, KSP oksp );   

};

class PetscCoupler {
  public : 
      static bool initialized;  
      static std::shared_ptr<PetscSS> static_petcs_ss_ptr; 
    
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

}

#endif 
#endif

