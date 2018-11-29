
#ifdef WITH_PETSCUI
#include "PetscInterface.h"

namespace SolverSelecter
{

PetscUI::PetscUI() : UserInterface<KSP,Vec>()
{
    srand(time(NULL));
    KSPCreate(PETSC_COMM_WORLD,&_ksp);
    internal_prefix = "internal_";
    
    AddParameter("mlinterface", "Waffles", "Choose a machine learning interface" );
    AddParameter("waffles.algorithm" , "RandomForest", "Choose the algorithm for Waffles");
    AddParameter("waffles.serialized", "", "Set the filename for waffles to read serialized model from");
    AddParameter("sqlite3.database", "./new_test_file.db", "Set the filename for sqlite3 database");
    AddParameter("C50.filestem", "./model", "set the filestem for the C50 serialized model");
    AddParameter("C50.database", "false", "build the model from the database at runtime (sloooow)");
    AddParameter("C50.trials", "-1", "Number of trials to use if building a C50 model");
    AddParameter("petsc.internalSample", "10", "Number of internal samples to use" );
    AddParameter("petsc.edgeSample", "10", "Number of edgeSamples to use" );
    AddParameter("petsc.matvec", "true", "Use Matvecs to extract samples" );
    AddParameter("power.interval","10000", "Interval in ns to sample power");

}

ErrorFlag PetscUI::Initialize() {
  SetTestingSpace();
}

ErrorFlag
PetscUI::SetTestingSpace()
{
    int durr = GetParameterAsInt("power.interval");
    testing_space.reset( new DefaultPetscTestingSpace(std::chrono::nanoseconds(durr) ));
    return error_flag;
}

ErrorFlag
PetscUI::GetDataBaseImplimentation(std::shared_ptr< DatabaseInterface > &database)
{
  
#ifdef WITH_SQLITE3
        database.reset( new SqlDatabase() );
        database->SetParameter("Database", GetParameter("sqlite3.database") );
#else
        std::cout << "PetscUI uses the Sqlite3 support by default. Either recompile with sqlite3" 
                  << "enabled, or override the PetscUI::GetDataBaseImplimentation function to set "
                  << "the database implementation" << std::endl;
         std::abort();
#endif
       
    return error_flag;
}

ErrorFlag
PetscUI::GetMachineLearningModel(std::shared_ptr< MachineLearningInterface > &machine)
{

   if ( GetParameter("mlinterface") == "Waffles") {
#ifdef WITH_WAFFLES   
    machine.reset( new WafflesInterface() ) ;
    machine->SetParameter("algorithm", GetParameter("waffles.algorithm") );
    machine->SetParameter("serialized", GetParameter("waffles.serialized") );
#else
    std::cout << "Compiled without support for waffles. Either recompile with waffles support or \
                  override the PetscUI::GetMachineLearningModel function to set the required \
                  machine learning interface." << std::endl;
    std::abort(); 
#endif    
  } else {
    machine.reset( new C50Interface() );
    machine->SetParameter("filestem", GetParameter("C50.filestem"));
    machine->SetParameter("trials", GetParameter("C50.trials"));
    machine->SetParameter("database", GetParameter("C50.database"));
  }
  return error_flag;
}

ErrorFlag
PetscUI::ChangeSolver( KSP &A, Vec &x, Vec &b , bool &change_it )
{
    PetscInt its = A->its;
    if ( its == 0 )
        change_it = true;
    else
        change_it = false;
    return 0;
}

ErrorFlag
PetscUI::SolveSystem( KSP &ksp, Vec &x, Vec &b, Solver &solver, bool &success)
{
    
    if (use_internal_ksp) 
    {
        PetscCopyFunction(_ksp, ksp );
        SetSolver(_ksp, solver );
        Mat AA,PP;
        KSPGetOperators(ksp, &AA, &PP );
        KSPSetOperators(_ksp, AA,PP);
        KSPSolve(_ksp, b, x );
        ksp->its = _ksp->its;
        ksp->totalits += ksp->its;
        ksp->reason = _ksp->reason;
    }
    else
    {
        SetSolver(ksp, solver);
        KSPSolve(ksp,b,x);
    }    
    success = (ksp->reason > 0 ) ;
    return 0;

}

int PetscUI::GetNNZ(KSP &ksp) {

      Mat A,P;
      MatInfo info;
      KSPGetOperators(ksp, &A, &P );
      MatGetInfo(A, MAT_GLOBAL_SUM, &info);

      PetscLogDouble nnz = info.nz_used;
      return (int) nnz;
}


ErrorFlag
PetscUI::SetSolver( KSP ksp, Solver solver )
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

    return error_flag;
}

ErrorFlag
PetscUI::ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap)
{
      
  testing_space->extract_features( ksp, fmap, GetParameterAsInt("petsc.edgeSample"), GetParameterAsInt("petsc.internalSample"), GetParameterAsBool("petsc.matvec") );
    return error_flag;
}

ErrorFlag
PetscUI::StartMeasurements(KSP &A,  Vec &x,  Vec &b )
{
    testing_space->start_measurements( A, x, b );
    return error_flag;
}

ErrorFlag
PetscUI::StopMeasurements( KSP &A, Vec &x, Vec &b, std::map<std::string, double> &mmap )
{
    testing_space->stop_measurements( A, x, b, mmap );
    return error_flag;
}

ErrorFlag
PetscUI::ClassificationMeasurements( std::map< std::string, double > &cmap )
{

    testing_space->classify_measurements( cmap );
    return error_flag;
}

ErrorFlag
PetscUI::InitMatrix( std::string filename, std::unique_ptr<KSP> &A ) 
{

    //Init matrix is called during database building when a new matrix is going
    //to be tested. So, we don't need to use the internal ksp. The external one 
    //is just fine. 
    use_internal_ksp=false;

    A.reset( new KSP() );
    KSPCreate(PETSC_COMM_WORLD,A.get());

    
    Mat mat_op;
    PetscViewer fd;

    MatCreate(PETSC_COMM_WORLD, &mat_op);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &fd);
    MatSetFromOptions(mat_op);
    MatLoad(mat_op,fd);
    
    PetscViewerDestroy(&fd);
    KSPSetOperators(*A, mat_op, mat_op);

    return error_flag;

}

ErrorFlag
PetscUI::InitVector( const KSP &A, std::unique_ptr<Vec> &x ) 
{

    x.reset( new Vec() );
    Mat mat_op, mat_pre;
    KSPGetOperators(A,&mat_op,&mat_pre);
    MatCreateVecs(mat_op, x.get(), NULL );
    VecAssemblyBegin(*x);
    VecAssemblyEnd(*x);

    return error_flag;
}

ErrorFlag
PetscUI::CopyVector( const Vec &x,
                     std::unique_ptr<Vec> &xx )
{
    xx.reset( new Vec() );
    VecDuplicate( x, xx.get());
    VecCopy(x,*xx);
    return error_flag;
}



ErrorFlag
PetscUI::FreeVector( std::unique_ptr<Vec> &x) 
{
    VecDestroy(x.get());
    return error_flag;
}

ErrorFlag
PetscUI::FreeMatrix( std::unique_ptr<KSP> &A ) 
{
    KSPDestroy(A.get());
    return error_flag;
}

ErrorFlag
PetscUI::SetVector( std::unique_ptr<Vec> &x, std::string type_ ) 
{

    if ( type_ == "ones" )
        VecSet(*x,1.0);
    else
        VecSetRandom(*x, NULL);

    VecAssemblyBegin(*x);
    VecAssemblyEnd(*x);
    return error_flag;
}

ErrorFlag
PetscUI::GetDefaultSolver( Solver &solver ) 
{
    solver.SetSolverName("gmres", "bjacobi");
    return error_flag;
}

PetscUI::~PetscUI()
{
    KSPDestroy(&_ksp);
}

void
PetscUI::PetscCopyFunction( KSP nksp, KSP oksp )
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
    nksp->pc_side             =   oksp->pc_side            ;
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
    nksp->pc_side_set         =   oksp->pc_side_set        ;
    nksp->normtype_set        =   oksp->normtype_set       ;
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
  PetscSS *ss = (PetscSS* )ksp->data;

    PetscFunctionBegin;
    Vec x,b;
    KSPGetSolution(ksp, &x );
    KSPGetRhs(ksp, &b );
    ss->Solve(ksp, x, b );
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

    PetscSS *ss = (PetscSS*) ksp->data;
    PetscUI *ui = (PetscUI*) ss->interface.get();
    ss->PrintSolver(ss->solver);
    KSPView( ui->_ksp, viewer);
    PetscFunctionReturn(0);
}

PetscErrorCode
PetscCoupler::KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
    PetscFunctionBegin;
    PetscSS *ss = (PetscSS*) ksp->data ;
    PetscUI *ui = (PetscUI*) ss->interface.get();
    bool flag = false;

    PetscOptionsHead(PetscOptionsObject,"KSP SS options");

    PetscBool flg;
    
    std::map<std::string, std::string> parametersMap;
    for ( auto &it: ui->GetParameters() ) {
        
        char inname[256];
        PetscStrcpy(inname, it.second.c_str());
        std::string name = "-ksp_ss_" + it.first;
        std::string helpp = ui->GetParameterHelp(it.first); 
        PetscOptionsString(name.c_str(), helpp.c_str(), NULL, inname,inname,256,&flg);
        if(flg) {
          parametersMap.insert(std::make_pair(it.first, inname));
        }
    } 
    ss->Initialize(parametersMap);
    PetscOptionsTail();
    PetscFunctionReturn(0);
}

// Decalare the static shared_ptr in the coupler
bool PetscCoupler::initialized = false; 
std::shared_ptr<PetscSS> PetscCoupler::static_petcs_ss_ptr = nullptr;

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

    KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,4);
    KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_RIGHT,3);
    KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_SYMMETRIC,2);
    KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_RIGHT,1);
    KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_LEFT,1);
   
    // We keep a copy of the shared pointer in a static member of the PetscCoupler class. That
    // way we do not have to rebuild the model every time the KSP is initialized.  
    if (PetscCoupler::static_petcs_ss_ptr == nullptr)
        PetscCoupler::static_petcs_ss_ptr.reset( new PetscSS( std::make_shared< PetscUI >() ));
    ksp->data = (void*) static_petcs_ss_ptr.get(); 
    
    // The solver selector should not use a preconditioner externally (although, we could actually 
    // set it up so that the preconditioner is actually a feature extraction?) 

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE); 
    PetscFunctionReturn(0);
}

}

#endif
