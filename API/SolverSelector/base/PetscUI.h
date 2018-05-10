#ifndef PETSCUI_HEADER
#define PETSCUI_HEADER

/** \file PetscSolverSelecter.h
 * \brief In this example, we precompile a solver selecter based on Petsc. This creates a library
 * that can be used in other example. This is different from the dummy example, where the entire
 * Solver selecter library is built at compilation time of the example.
 **/


//Include the petsc files
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "SolverSelecter.h"
#include "PetscTestingSpace.h"

/* First we need to define the Petsc User Interface */
class PetscUI : public _SS_UserInterface<KSP,Vec>
{
public:

    KSP _ksp; /**< linear solver */
    std::string binary_file;

    std::shared_ptr< PetscTestingSpace > testing_space;
   
    PetscUI() : _SS_UserInterface<KSP,Vec>() 
    {
        srand(time(NULL));
        KSPCreate(PETSC_COMM_WORLD,&_ksp);
        SetTestingSpace("default");
    }


    /* The interface functions that must be overriden. Some are optional, but, we impliment them all */
    _SS_ErrorFlag SolveSystem( KSP &ksp, Vec &x, Vec &b, _SS_Solver &solver) override
    {

        /* In this function, we copy over the information to a new, local ksp. This occurs
        * because we cannot change the input ksp type at this point. ( This function is
        * called from inside the KSPSolve function of the solver-selecter.) */

        std::string solvern, precond;
        std::set< std::string >  keys;
        solver.GetSolverInfo( solvern, precond, keys );

        /* Set some defualt parameters, just in case */
        if ( solvern == "NONE" )
            solvern = "gmres";
        if ( precond == "NONE" )
            precond = "bjacobi";
   
        /* Set the solver and tolerances */
        KSPSetType(_ksp, solvern.c_str() );
        KSPSetTolerances(_ksp,1e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

        /* Set the Preconditioner options */
        PC pc;
        KSPGetPC(_ksp,&pc);
        PCSetType(pc,precond.c_str());

        /* Set all the petsc options */
        std::string pvalue;
        if ( keys.size() > 0 )
        {
            for ( auto it : keys )
            {
                solver.GetParameter(it, pvalue);
                PetscOptionsSetValue(NULL, it.c_str(), pvalue.c_str());
            }
        }
        /* Set the operators */
        Mat AA,PP;
        KSPGetOperators(ksp, &AA, &PP);
        KSPSetOperators(_ksp, AA, PP );

        /* Solve the system */
        KSPSolve(_ksp, b, x );

        /* Need to set the convergence critera for the outer SS KSP type */
        ksp->reason   = _ksp->reason;
        ksp->its      = _ksp->its;
        ksp->totalits = _ksp->totalits;

        return 0;
    }

    _SS_ErrorFlag SetTestingSpace( std::string space ) {
          
         /* Add in testing spaces here */ 
         if ( space == "default") {
             testing_space.reset( new DummyTestingSpace() );
         else {
             testing_space.reset( new DummyTestingSpace() ); 
         }
    }

    _SS_ErrorFlag ExtractFeatures( KSP &ksp, std::map<std::string, double> &fmap) override
    {
        testing_space.extract_features( ksp, fmap );
        return _SS_error_flag;
    }
 
    _SS_ErrorFlag StartMeasurements( const KSP &A, const Vec &x, const Vec &b ) override
    {
        testing_space.start_measurements( A, x, b );
        return _SS_error_flag;
    }

    _SS_ErrorFlag StopMeasurements( const  KSP &A, const Vec &x, const Vec &b, std::map<std::string, double> &mmap ) override {
        testing_space.stop_measurements( A, x, b, mmap );
        return _SS_error_flag;
    }

    _SS_ErrorFlag ClassificationMeasurements( std::map< std::string, double > &cmap ) override {        

        testing_space.classify_measurements( cmap ); 
        return _SS_error_flag;
    }

    _SS_ErrorFlag GetMachineLearningModel(std::shared_ptr< _SS_MachineLearning > &machine)
    {
        testing_space.set_machine_learning_model( machine );
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

        PetscInt n,m;
        Mat mat_op, mat_pre;
        KSPGetOperators(A,&mat_op,&mat_pre);
        MatGetSize(mat_op, &m, &n);

        x.reset( new Vec() );
        VecCreate(MPI_COMM_WORLD, x.get() );
        VecSetSizes(*x, PETSC_DECIDE, m );
        VecSetFromOptions(*x);
        VecAssemblyBegin(*x);
        VecAssemblyEnd(*x);
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
        /* check the vector size and cols are the same */
        PetscInt size;
        VecGetSize(*x, &size);

        //*x = (PetscScalar) 0.0;
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



    _SS_ErrorFlag Finalize()
    {
        KSPDestroy(&_ksp);
        return _SS_error_flag;
    } 
    
};

class PetscSS : public _SS_SolverSelecter<KSP,Vec>
{
public:
    PetscSS( std::shared_ptr< PetscUI > _interface ) : _SS_SolverSelecter<KSP,Vec>( _interface )
    {
    }

    /* To use the solver selecter we need to wrap it up inside a KSPType. The reaminder of the class
     * is Petsc boiler plate to do just that. Simply call */
    _SS_ErrorFlag WrapAsKSP( KSP &ksp )
    {
        /* Register the Petsc Solver Selecter */
        KSPRegister("KSPSS", KSPCreate_SS);
        KSPSetType(ksp, "KSPSS");

        /* Get the preconditioner */
        PC pcc;
        KSPGetPC(ksp,&pcc);
        PCSetType(pcc,PCNONE);

        /* Set up the KSP solver */
        KSPSetFromOptions(ksp); /* Not impliented yet */
        KSPSetSkipPCSetFromOptions(ksp,PETSC_TRUE);
        KSPSetUp(ksp); /* This builds the data structures */

        return 0;
    }

    static PetscErrorCode KSPSetUp_SS(KSP ksp)
    {
        PetscFunctionBegin;
        ksp->data = (void*) this ;        
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
        PetscFunctionReturn(0);
    }

    static PetscErrorCode KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp)
    {
        PetscFunctionBegin;
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
        PetscFunctionReturn(0);
    }

};

#endif
