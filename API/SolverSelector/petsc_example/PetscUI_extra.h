#ifndef HEADER_PETSCUI_EXTRA
#define HEADER_PETSCUI_EXTRA

//Include the petsc files
#include <cstring>
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "typedefs.h"

class PetscTestingSpace {
  public:
  PetscTestingSpace( ) { }
  virtual void extract_features( KSP &ksp, std::map<std::string, double> &fmap ) = 0;
  virtual void start_measurements( KSP &ksp, Vec &x, Vec &b ) = 0;
  virtual void stop_measurements(KSP &ksp, Vec &x, Vec &b , std::map<std::string, double > &mmap) = 0;
  virtual void classify_measurements( std::map<std::string, double> &cmap ) = 0;
  virtual void set_machine_learning_model ( std::shared_ptr<_SS_MachineLearning> &machine ) = 0;
};

void PetscCopyFunction( KSP nksp, KSP oksp ) {


  nksp->dm                  =   oksp->dm                 ; 
  nksp-> dmAuto             =   oksp->dmAuto             ;
  nksp->dmActive            =   oksp->dmActive           ;
  nksp->max_it              =   oksp->max_it             ;
  nksp->guess               =   oksp->guess              ;
  nksp->guess_zero         =   oksp->guess_zero        ;
  nksp->calc_sings         =   oksp->calc_sings        ;
  nksp->calc_ritz          =   oksp->calc_ritz         ;
  nksp->guess_knoll         =   oksp->guess_knoll        ;
  nksp->pc_side             =   oksp->pc_side            ;
  nksp->rtol               =   oksp->rtol              ;
  nksp->abstol             =   oksp->abstol            ;
  nksp->ttol               =   oksp->ttol              ;
  nksp->divtol              =   oksp->divtol             ;
  nksp->rnorm0              =   oksp->rnorm0             ;
  nksp->rnorm               =   oksp->rnorm              ;
  nksp->reason              =   oksp->reason             ;
  nksp->errorifnotconverged =   oksp->errorifnotconverged;
  //nksp->vec_sol            =   oksp->vec_sol           ;
  //nksp->vec_rhs             =   oksp->vec_rhs            ;
  nksp->res_hist            =   oksp->res_hist           ;
  nksp->res_hist_alloc      =   oksp->res_hist_alloc     ;
  nksp->res_hist_len        =   oksp->res_hist_len       ;
  nksp->res_hist_max        =   oksp->res_hist_max       ;
  nksp->res_hist_reset      =   oksp->res_hist_reset     ;
  nksp->chknorm             =   oksp->chknorm            ;
  nksp->lagnorm             =   oksp->lagnorm            ;  
  
  
  nksp->numbermonitors      =   oksp->numbermonitors     ;
  nksp->converged           =   oksp->converged          ;
  nksp->convergeddestroy    =   oksp->convergeddestroy   ;
  nksp->user                =   oksp->user               ;  
  
  //nksp->setupstage          =   oksp->setupstage         ;
  //nksp->its                 =   oksp->its                ;
  //nksp->totalits            =   oksp->totalits           ;
  nksp->transpose_solve     =   oksp->transpose_solve    ;
  //nksp->normtype            =   oksp->normtype           ;
  nksp->pc_side_set         =   oksp->pc_side_set        ;
  nksp->normtype_set        =   oksp->normtype_set       ;
  nksp->dscale              =   oksp->dscale             ;
  nksp->dscalefix           =   oksp->dscalefix          ;
  nksp->dscalefix2          =   oksp->dscalefix2         ;
  nksp->diagonal            =   oksp->diagonal           ;
  nksp->truediagonal        =   oksp->truediagonal       ;
  nksp->eigviewer           =   oksp->eigviewer          ; 
  
  //Arrays that need to be copied 
  //std::memcpy(nksp->normsupporttable, oksp->normsupporttable , sizeof(nksp->normsupporttable) );
  std::memcpy(nksp->monitor         , oksp->monitor          , sizeof(nksp->monitor         ) );
  std::memcpy(nksp->monitordestroy  , oksp->monitordestroy   , sizeof(nksp->monitordestroy  ) );
  std::memcpy(nksp->monitorcontext  , oksp->monitorcontext   , sizeof(nksp->monitorcontext  ) );         
   
  // (new version ) nksp->setupnewmatrix      =   oksp->setupnewmatrix     ;
  // (new version ) nksp->skippcsetfromoptions=   oksp->skippcsetfromoption;
}                           
                  
#endif


