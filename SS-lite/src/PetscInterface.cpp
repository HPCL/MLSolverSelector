



#include "PetscInterface.h"
#include "PetscFeatureSetSelector.h"


PetscInterface::PetscInterface(){

    //Create the internal ksp used for the solves (might be redundent ? ) 
    KSPCreate(PETSC_COMM_WORLD,&_ksp);
    internal_prefix = "internal_";
};

int PetscInterface::Initialize(std::string filestem_, std::string featureSet_, int edge, int interior, bool matvec) {

  if ( initialized ) return 1;  
  
  filestem = filestem_;
    featureSet = featureSet_;
    edgeSamples = edge;
    interiorSamples = interior;
    useMatVecs = matvec;
    
    //Create the C50 Interface for the given filestem
    
    c50Interface.reset( new C50Interface(filestem) ); 
    

    if ( !PetscFeatureSetSelector::featureSetExists(featureSet) ) {
        std::cout << " A feature set with the name " << featureSet << "does not exist \n" 
                  << " You need to go into the feature sets directory, create a new " 
                  << featureSet << ".fs file with the correct features in it, then " 
                  << "go back one directory and run the ./generate... shell script, then"  
                  << " recompile (make clean && make). This is kind of annoying, but it ensures that the user"  
                  << " at least thinks about defining a feature set that matches the features" 
                  << " in the C50 model \n";
        std::abort();
    } 
    
    initialized = true;
    return 1;
}

int 
PetscInterface::PrintSolver() {
  currentSolver.Print();
  return 1;
}

int PetscInterface::ClassifyAndSolve(KSP &ksp) {
      
    bool success;
    std::map<std::string, double> featuresMap;
    
    ExtractFeatures(ksp, featuresMap);
    
    std::cout << " SolverSelector Extracted the following Features from the Matrix : \n";
    for ( auto it : featuresMap ) {
        std::cout << " \t " << it.first << " : " << it.second << "\n";
    }

    
    Predict(featuresMap, currentSolver);
    
    std::cout << " Solver Selector Selected The following Solver For this Matrix \n";
    currentSolver.Print();

    
    SolveSystem(ksp, currentSolver, success);
    
    PetscInt its;
    KSPGetIterationNumber(ksp,&its);  
    if ( success > 0 )  { 
      std::cout << " Success -- Solve converged reason (see petsc KSPConverged Reason : " << success << "\n" ; 
      std::cout << " Iteration Count : " << its ;                                                                                               
    }
    else  {
          std::cout << " The chosen solver didn't work: " << success << "\n";
          std::cout << " Using the default instead \n";
          GetDefaultSolver(currentSolver);
          SolveSystem(ksp, currentSolver, success); 
          if ( !success ) {
            printf("The default solver failed as well. ¯\\_(ツ)_/¯ \n");
          }
        }
    return 1;  
}
              
int PetscInterface::Predict(std::map<std::string, double> &features_map, Solver &solver) {
    solver.Clear();
    solver = c50Interface->Predict(features_map);
    return 1;
}

int
PetscInterface::SolveSystem( KSP &ksp, Solver &solver, bool &success)
{

  

    PetscCopyFunction(_ksp, ksp );
    SetSolver(_ksp, solver );
    Mat AA,PP;
    KSPGetOperators(ksp, &AA, &PP );
    KSPSetOperators(_ksp, AA,PP);
    Vec x,b;
    KSPGetSolution(ksp, &x );
    KSPGetRhs(ksp, &b );
    KSPSolve(_ksp, b, x );
    ksp->its = _ksp->its;
    ksp->totalits += ksp->its;
    ksp->reason = _ksp->reason;
    success = (ksp->reason > 0 ) ;
    return 0;

}

int
PetscInterface::SetSolver( KSP ksp, Solver solver )
{


    std::string solvern, precond;
    std::set< std::string >  keys;
    solver.GetSolverInfo( solvern, precond, keys );
    if (precond == "NONE") precond = "none";

    /* Set the solver and Preconditioner */
    PC pc;
    KSPGetPC(ksp,&pc);
    KSPSetType(ksp, solvern.c_str() );
    PCSetType(pc,precond.c_str());

    std::string pvalue;
    if ( keys.size() > 0 )
    {
        for ( auto it : keys )
        {
            solver.GetParameter(it, pvalue);
            std::string pname = "-" + internal_prefix + it;
            PetscOptionsSetValue(NULL, pname.c_str(), pvalue.c_str());
        }
    }

    KSPSetOptionsPrefix(ksp, internal_prefix.c_str());
    KSPSetFromOptions(ksp);

    return 0;
}

int
PetscInterface::ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap)
{
   Mat AA,PP;                                                                                                
   KSPGetOperators(ksp, &AA, &PP );                                                                          
   PetscFeatureSetSelector::ExtractJacobianFeatures(featureSet,AA,edgeSamples,interiorSamples,fmap,useMatVecs);   
   
   
   
   return 0;
}

int
PetscInterface::GetDefaultSolver( Solver &solver ) 
{
    solver.SetSolverName("gmres", "bjacobi");
    return 0;
}

PetscInterface::~PetscInterface()
{
    KSPDestroy(&_ksp);
}

void
PetscInterface::PetscCopyFunction( KSP nksp, KSP oksp )
{

    nksp->dm                  =   oksp->dm                 ;
    nksp-> dmAuto             =   oksp->dmAuto             ;
    nksp->dmActive            =   oksp->dmActive           ;
    nksp->max_it              =   oksp->max_it             ;
    nksp->guess               =   oksp->guess              ;
    nksp->guess_zero          =   oksp->guess_zero        ;
    nksp->calc_sings          =   oksp->calc_sings        ;
    nksp->calc_ritz           =   oksp->calc_ritz         ;
    nksp->guess_knoll         =   oksp->guess_knoll        ;
    //nksp->pc_side             =   oksp->pc_side            ;
    nksp->rtol                =   oksp->rtol              ;
    nksp->abstol              =   oksp->abstol            ;
    nksp->ttol                =   oksp->ttol              ;
    nksp->divtol              =   oksp->divtol             ;
    nksp->rnorm0              =   oksp->rnorm0             ;
    nksp->rnorm               =   oksp->rnorm              ;
    nksp->reason              =   oksp->reason             ;
    nksp->errorifnotconverged =   oksp->errorifnotconverged;
    nksp->res_hist            =   oksp->res_hist           ;
    nksp->res_hist_alloc      =   oksp->res_hist_alloc     ;
    nksp->res_hist_len        =   oksp->res_hist_len       ;
    nksp->res_hist_max        =   oksp->res_hist_max       ;
    nksp->res_hist_reset      =   oksp->res_hist_reset     ;
    nksp->chknorm             =   oksp->chknorm            ;
    nksp->lagnorm             =   oksp->lagnorm            ;
    nksp->cnvP                =   oksp->cnvP               ;

    nksp->numbermonitors      =   oksp->numbermonitors     ;
    nksp->converged           =   oksp->converged          ;
    nksp->convergeddestroy    =   oksp->convergeddestroy   ;
    nksp->user                =   oksp->user               ;
    nksp->transpose_solve     =   oksp->transpose_solve    ;
    //nksp->pc_side_set         =   oksp->pc_side_set        ;
    //nksp->normtype_set        =   oksp->normtype_set       ;
    nksp->dscale              =   oksp->dscale             ;
    nksp->dscalefix           =   oksp->dscalefix          ;
    nksp->dscalefix2          =   oksp->dscalefix2         ;
    nksp->diagonal            =   oksp->diagonal           ;
    nksp->truediagonal        =   oksp->truediagonal       ;
    nksp->eigviewer           =   oksp->eigviewer          ;

    //Arrays that need to be copied
    std::memcpy(nksp->monitor         , oksp->monitor          , sizeof(nksp->monitor         ) );
    std::memcpy(nksp->monitordestroy  , oksp->monitordestroy   , sizeof(nksp->monitordestroy  ) );
    std::memcpy(nksp->monitorcontext  , oksp->monitorcontext   , sizeof(nksp->monitorcontext  ) );

}

PetscErrorCode
PetscCoupler::KSPSetUp_SS(KSP ksp)
{
    PetscFunctionBegin;
    PetscFunctionReturn(0);
}

PetscErrorCode
PetscCoupler::KSPSolve_SS(KSP ksp)
{
  
    PetscInterface *ss = (PetscInterface* )ksp->data;
    PetscFunctionBegin;
    ss->ClassifyAndSolve(ksp);
    PetscFunctionReturn(0);
}

PetscErrorCode
PetscCoupler::KSPDestroy_SS(KSP ksp)
{
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}

PetscErrorCode
PetscCoupler::KSPView_SS(KSP ksp, PetscViewer viewer)
{
    PetscFunctionBegin;
    PetscInterface *ss = (PetscInterface*) ksp->data;
    KSPView( ss->_ksp, viewer);
    PetscFunctionReturn(0);
}

PetscErrorCode
PetscCoupler::KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
    PetscFunctionBegin;
    PetscInterface *ss = (PetscInterface*) ksp->data ;
    bool flag = false;

    PetscOptionsHead(PetscOptionsObject,"KSP SS-lite options");

    PetscBool flg;
    
   
    PetscBool matvecs = PETSC_TRUE;
    PetscOptionsBool("-ksp_ss_matvecs", "Should we use Matvecs to extract columns", NULL, matvecs, &matvecs, &flg);
    PetscInt edge = 10;
    PetscOptionsInt("-ksp_ss_edge", "Number of edge points to sample", NULL, edge, &edge, &flg);
    PetscInt inter = 10;
    PetscOptionsInt("-ksp_ss_interior", "Number of interior points to sample", NULL, inter, &inter, &flg);
    char filestem[256] = ""; 
    PetscOptionsString("-ksp_ss_filestem", "C5.0 file stem", NULL, filestem, filestem,256,&flg);
    char fset[256] = ""; 
    PetscOptionsString("-ksp_ss_featureset", "Name of feature set", NULL, fset, fset,256,&flg);
    
    ss->Initialize( std::string(filestem), std::string(fset), edge, inter, matvecs);
    PetscOptionsTail();
    PetscFunctionReturn(0);
}

// Decalare the static shared_ptr in the coupler
bool PetscCoupler::initialized = false; 
std::shared_ptr<PetscInterface> PetscCoupler::static_petcs_interface_ptr = nullptr;

PetscErrorCode
PetscCoupler::KSPCreate_SS(KSP ksp)
{
    PetscFunctionBegin;
    ksp->ops->setup          = KSPSetUp_SS;
    ksp->ops->solve          = KSPSolve_SS;
    ksp->ops->destroy        = KSPDestroy_SS;
    ksp->ops->view           = KSPView_SS;
    ksp->ops->setfromoptions = KSPSetFromOptions_SS;
    ksp->ops->buildsolution  = KSPBuildSolutionDefault;
    ksp->ops->buildresidual  = KSPBuildResidualDefault;
 
    // We keep a copy of the shared pointer in a static member of the PetscCoupler class. That
    // way we do not have to reload the model every time the KSP is initialized.  
    if (PetscCoupler::static_petcs_interface_ptr == nullptr)
        PetscCoupler::static_petcs_interface_ptr.reset( new PetscInterface());
    
    ksp->data = (void*) static_petcs_interface_ptr.get(); 
    
    // The solver selector should not use a preconditioner externally (although, we could actually 
    // set it up so that the preconditioner is actually a feature extraction?) 

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE); 
    PetscFunctionReturn(0);
}

