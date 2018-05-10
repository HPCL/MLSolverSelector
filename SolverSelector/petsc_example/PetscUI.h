#ifndef PETSCUI_HEADER
#define PETSCUI_HEADER

/** \file PetscSolverSelecter.h
 * \brief In this example, we precompile a solver selecter based on Petsc. This creates a library
 * that can be used in other example. This is different from the dummy example, where the entire
 * Solver selecter library is built at compilation time of the example.
 **/

#include "SolverSelecter.h"
#include "PetscUI_extra.h"
#include "testing_spaces/DummyTestingSpace.h"

/* First we need to define the Petsc User Interface */
class PetscUI : public _SS_UserInterface<KSP,Vec>
{
public:

    KSP _ksp; /**< linear solver */
    std::shared_ptr< PetscTestingSpace > testing_space;

    bool pre_extract;
    std::map< std::string, double> pre_fmap; /* Pre extracted features calculated during setup */
    std::string internal_prefix;
    PetscUI() : _SS_UserInterface<KSP,Vec>() 
    {
        srand(time(NULL));
        KSPCreate(PETSC_COMM_WORLD,&_ksp);
        SetTestingSpace("default");
        input_file = "../inputs/petsc.input" ;
        database_name = "../databases/petsc.database";
        internal_prefix = "internal_";
    }


    /* The interface functions that must be overriden. Some are optional, but, we impliment them all */
    _SS_ErrorFlag SolveSystem( KSP &ksp, Vec &x, Vec &b, _SS_Solver &solver) override
    {

        if ( ! pre_extract ) {
            SetInternalKSP( _ksp, ksp, solver );
        }
        KSPSolve(_ksp, b, x );       
        ksp->its = _ksp->its;
        ksp->totalits += ksp->its;
        ksp->reason = _ksp->reason;
        
        return 0;
    }

    _SS_ErrorFlag SetInternalKSP( KSP ksp, _SS_Solver solver ) {
      
        SetInternalKSP( _ksp, ksp, solver );
        return _SS_error_flag;
    }

    _SS_ErrorFlag SetInternalKSP( KSP nksp, KSP oksp, _SS_Solver solver ) {
        
        /* Reset the internal KSP Solver to be a soft clone of the external solver */ 
        PetscCopyFunction(nksp, oksp); 
        
        std::string solvern, precond;
        std::set< std::string >  keys;
        solver.GetSolverInfo( solvern, precond, keys ); 
        if (precond == "NONE") precond = "none";
        
        /* Set the solver and Preconditioner */
        PC pc;
        KSPGetPC(nksp,&pc);
        KSPSetType(nksp, solvern.c_str() );
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
        
        KSPSetOptionsPrefix(_ksp, internal_prefix.c_str());
        KSPSetFromOptions(_ksp);
        
        Mat AA,PP;
        KSPGetOperators(oksp, &AA, &PP );
        KSPSetOperators(nksp, AA,PP);
    }

 
    _SS_ErrorFlag SetTestingSpace( std::string space ) {
          
         /* Add in testing spaces here */ 
         if ( space == "default") testing_space.reset( new DummyTestingSpace() );
         else testing_space.reset( new DummyTestingSpace() ); 
         
    }

    _SS_ErrorFlag PreExtractFeatures( KSP &ksp ) {
        
        pre_fmap.clear();
        testing_space->extract_features( ksp, pre_fmap );
        return _SS_error_flag;
    }

    _SS_ErrorFlag ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap) override
    {
        std::cout << "PRINTING" << std::endl;
        
        Mat AA,PP;
         
        KSPGetOperators(ksp, &AA, &PP );
        
        if ( pre_extract && pre_fmap.size() > 0 ) fmap = pre_fmap;  
        else testing_space->extract_features( ksp, fmap );
        return _SS_error_flag;
    }
 
    _SS_ErrorFlag StartMeasurements(KSP &A,  Vec &x,  Vec &b ) override
    {
        testing_space->start_measurements( A, x, b );
        return _SS_error_flag;
    }

    _SS_ErrorFlag StopMeasurements( KSP &A, Vec &x, Vec &b, std::map<std::string, double> &mmap ) override {
        testing_space->stop_measurements( A, x, b, mmap );
        return _SS_error_flag;
    }

    _SS_ErrorFlag ClassificationMeasurements( std::map< std::string, double > &cmap ) override {        

        testing_space->classify_measurements( cmap ); 
        return _SS_error_flag;
    }

    _SS_ErrorFlag GetMachineLearningModel( std::shared_ptr<_SS_MachineLearning> &machine)
    {
        testing_space->set_machine_learning_model( machine );
        return _SS_error_flag;
    }    

    _SS_ErrorFlag InitMatrix( std::string filename, std::unique_ptr<KSP> &A ) override
    {

        /* Init matrix is potentially called on multiple different matricies when
         * building the database, so, we need to delete the KSP solver, each time
         * a new matrix is used. This is never called during solver selection, only
         * during a run where we are building the database. */
        KSPDestroy( &_ksp);
        KSPCreate(PETSC_COMM_WORLD, &_ksp);

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
        //MatView(mat_op,PETSC_VIEWER_STDOUT_WORLD);

        return _SS_error_flag;

    }

    _SS_ErrorFlag InitVector( const KSP &A, std::unique_ptr<Vec> &x ) override
    {
 
        x.reset( new Vec() );
        Mat mat_op, mat_pre;
        KSPGetOperators(A,&mat_op,&mat_pre); 
        MatCreateVecs(mat_op, x.get(), NULL );
        VecAssemblyBegin(*x);
        VecAssemblyEnd(*x);
 
        return _SS_error_flag;
    }

    virtual _SS_ErrorFlag CopyVector( const Vec &x,
                                      std::unique_ptr<Vec> &xx )
    {
        xx.reset( new Vec() );
        VecDuplicate(x,xx.get());
        return _SS_error_flag;
    }



    _SS_ErrorFlag FreeVector( std::unique_ptr<Vec> &x) override
    {
        VecDestroy(x.get());
        return _SS_error_flag;
    }

    _SS_ErrorFlag FreeMatrix( std::unique_ptr<KSP> &A ) override
    {
        KSPDestroy(A.get());
        return _SS_error_flag;
    }

    _SS_ErrorFlag SetVector( std::unique_ptr<Vec> &x, std::string type_ ) override
    {

        if ( type_ == "ones" )
            VecSet(*x,1.0);
        else
            VecSetRandom(*x, NULL);

        VecAssemblyBegin(*x);
        VecAssemblyEnd(*x);
        return _SS_error_flag;
    }

    _SS_ErrorFlag GetDefaultSolver( _SS_Solver &solver ) override
    {
        solver.SetSolverName("gmres", "NONE");
        return _SS_error_flag;
    }


    virtual  ~PetscUI()
    {
        KSPDestroy(&_ksp);
    } 
    
};

typedef _SS_SolverSelecter<KSP,Vec> PetscSS;


static PetscErrorCode KSPSetUp_SS(KSP ksp)
{
    PetscFunctionBegin;
   
    /* Make sure the preconditioner is turned off */ 
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
  
    /* For PETSc, we do things a little different. Normally, you would 
     * just call the ss->Solve function, but for Petsc we need to do things
     * in parts. */
        
    PetscSS *ss = (PetscSS*) ksp->data;
    PetscUI *interface = (PetscUI*) ss->interface.get();
    _SS_Solver solver;
    if ( interface->pre_extract  ) {
        if ( interface->build_inline || interface->use_ml )
        {
            interface->GetDefaultSolver(solver) ;
        }
        else
        {
            std::map<std::string, double> features_map;
            ss->interface->ExtractFeatures( ksp, features_map);
            ss->machinemodel->Classify(features_map, solver );

            if (solver.solver == "NONE" )
                interface->GetDefaultSolver(solver);
        }
        interface->SetInternalKSP( ksp, solver );
        return _SS_error_flag;
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode KSPSolve_SS(KSP ksp)
{
    PetscSS *ss = ( PetscSS* )ksp->data;

    PetscFunctionBegin;
    Vec x,b;
    KSPGetSolution(ksp, &x );
    KSPGetRhs(ksp, &b );
    ss->Solve(ksp, x, b );
    PetscFunctionReturn(0);
}

static PetscErrorCode KSPDestroy_SS(KSP ksp)
{
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}

static PetscErrorCode KSPView_SS(KSP ksp, PetscViewer viewer)
{
    PetscFunctionBegin;

    PetscSS *ss = (PetscSS*) ksp->data;
    ss->PrintSolver(ss->solver);
    PetscFunctionReturn(0);
}

static PetscErrorCode KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
    PetscFunctionBegin;
 
    PetscSS *ss = (PetscSS*) ksp->data ;
    PetscUI *ui = (PetscUI*) ss->interface.get();
    bool flag = false;

    PetscOptionsHead(PetscOptionsObject,"KSP SS options");
    
    PetscBool temp;

    temp = (PetscBool) ui->pre_extract;
    PetscOptionsBool("-ksp_ss_pre_extract","Pre extract features",NULL,temp,&temp,NULL);
    ui->pre_extract = (bool) temp;

    temp = (PetscBool) ui->use_ml;
    PetscOptionsBool("-ksp_ss_use_ml","Use machine learning",NULL,temp,&temp,NULL);
    ui->use_ml = (bool) temp;
     
    temp = (PetscBool) ui->build_inline;
    PetscOptionsBool("-ksp_ss_build_inline","Build database inline",NULL,temp, &temp ,NULL);
    ui->build_inline = (bool) temp;
   
    PetscBool flg;
    char dbname[256];
    PetscStrcpy(dbname, ui->database_name.c_str());
    PetscOptionsString("-ksp_ss_database", " Name of the database", "", dbname,dbname,256,&flg);
    if (flg) ui->database_name = dbname;

    char inname[256];
    PetscStrcpy(dbname, ui->input_file.c_str());
    PetscOptionsString("-ksp_ss_inputfile", " Name of the inputfile", NULL, inname,inname,256,&flg);
    if (flg) ui->input_file = inname;
  
    PetscOptionsTail();

    PetscFunctionReturn(0);
}

static PetscErrorCode KSPCreate_SS(KSP ksp)
{
    PetscFunctionBegin;
    ksp->ops->setup          = KSPSetUp_SS;
    ksp->ops->solve          = KSPSolve_SS;
    ksp->ops->destroy        = KSPDestroy_SS;
    ksp->ops->view           = KSPView_SS;
    ksp->ops->setfromoptions = KSPSetFromOptions_SS;
    ksp->ops->buildsolution  = KSPBuildSolutionDefault;
    ksp->ops->buildresidual  = KSPBuildResidualDefault;
    
    PetscSS *ss = new PetscSS( std::make_shared< PetscUI >() );
    ksp->data = (void*) ss;
    
    PetscFunctionReturn(0);
}



#endif