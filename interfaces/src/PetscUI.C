
#ifdef WITH_PETSCUI

#include "PetscUI.h"

using namespace petscUI;

/** These are the specific implimentations used in this petsc solver. This
 * should make it easy to test a whole bunch of features and machine learning
 * impls if we want to. */

/** \file PetscUI.C
 * \brief Implimentations for the PetscUI 
 **/

PetscUI::PetscUI(std::string _db, std::string _in, bool _so, bool _do, bool _io) : _SS_UserInterface<KSP,Vec>()
{
    dname = _db;
    inputfile = _in;
    solver_on = _so;
    dump_on = _do;
    inline_on = _io;
    srand(time(NULL));
    KSPCreate(PETSC_COMM_WORLD,&_ksp);
}

/** Add Features */
_SS_ErrorFlag
PetscUI::AddFeaturesAndMeasurements( _SS_Features &features, _SS_Measurements &measure )
{
    /* Several features are included in the features directory. Custom features are simple
     * to impliment based on those examples */
    std::shared_ptr<_SS_Feature> inf( new Features::Norm(-1, "InfNorm") );
    std::shared_ptr<_SS_Feature> fro( new Features::Norm(0, "FroNorm") );
    std::shared_ptr<_SS_Feature> var( new Features::Variance(1, "RowVar") );
    features.AddFeature(inf);
    features.AddFeature(fro);
    features.AddFeature(var);

    /* Only one measurement is currently supplied in the measurements directory, although
     * custom measurements are easy to write. One example might be memory consumption or
     * even accuracy, iteration counts, resilliancy, anything really. */
    std::shared_ptr<_SS_Measurement> cm(  new Measurements::CPUTime(0.3) );
    measure.AddMeasurement("CPUTime" , cm);
    return _SS_error_flag;
}

_SS_ErrorFlag 
PetscUI::GetDataBaseImpl( std::shared_ptr< _SS_DataBaseBase > &database )
{
  database.reset( new Database::_SS_DataBaseSql( dname ) );
  
  return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::GetMLImpl( std::shared_ptr< _SS_MachineLearning > &machinemodel )
{
  machinemodel.reset( new MachineLearning::_ML_Waffles() );
  return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::SetSolverSelecterOptions( _SS_SolverSelecter<KSP,Vec> *ss ) 
{

    ss->Set("database_name",dname);
    ss->SetSolverSelection(solver_on);
    ss->SetMatrixDump(dump_on);
    ss->SetBuildDatabaseInline(inline_on, inputfile );
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::GetMatrixInfo( KSP &A, int &nrows, int &ncols, int &chunks,
                                    std::string &matrix_name, bool &mfree ) 
{

    /* Set the matrix size ( Global ) */
    Mat AA,PP;
    KSPGetOperators(A, &AA, &PP);
    MatGetSize(AA,&nrows,&ncols);

    chunks = nrows; /* get the sparcity pattern row by row. This is SUPER slow. but, easy on the memory */
    std::ostringstream oss;
    oss << " inline_" << rand();
    matrix_name = oss.str();
    mfree = false;
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::GetSparcity( KSP &ksp, int &chunk,
                                  std::vector< std::pair< unsigned int, unsigned int > > &sparcity,
                                  std::vector< double > &values
                                )
{
    Mat AA,PP;
    KSPGetOperators(ksp, &AA, &PP);

    int nnz = 0;
    const PetscScalar *vals;
    const PetscInt *cols;

    MatGetRow(AA, chunk, &nnz, &cols, &vals);    
    sparcity.reserve(nnz);
    values.reserve(nnz);
    for (int j = 0; j < nnz; j++ )
    {
        sparcity.push_back( std::make_pair( chunk, (unsigned int) cols[j] ) );
        values.push_back( (double) vals[j] );
    }
    MatRestoreRow(AA, chunk, &nnz, &cols, &vals);
    return _SS_error_flag ;
}

_SS_ErrorFlag
PetscUI::MatVecAndGetData( KSP &ksp, Vec &x, std::vector<int> &row, std::vector<double> &val)
{

    Mat A, P;
    Vec y;
    KSPGetOperators(ksp,&A,&P);
    VecDuplicate(x,&y);
    MatMult(A,x,y);
    VecGetValues( y, row.size(), &row[0],&val[0]);
    VecDestroy(&y);
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::SolveSystem( KSP &ksp, Vec &x, Vec &b, _SS_Solver &solver, std::map<std::string,double> &mstruct) 
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

    /* Need to set the parameters in the mstruct that must be set */
    for (auto &it : mstruct )
    {
        if ( it.first == "converged" )
            it.second = (double) _ksp->reason;
        else if ( it.first == "another_parameter")
            it.second = 1.3333;
        else
            std::cout << " Woops, forgot to set " << it.first << " in the SolveSystem mstruct \n";
    }
    /* Need to set the convergence critera for the outer SS KSP type */
    ksp->reason   = _ksp->reason;
    ksp->its      = _ksp->its;
    ksp->totalits = _ksp->totalits;

    return 0;
}


_SS_ErrorFlag
PetscUI::InitMatrix( std::string filename, std::unique_ptr<KSP> &A ) 
{

    /* Init matrix is potentially called on multiple different matricies when
     * building the database, so, we need to delete the KSP solver, each time
     * a new matrix is used */
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

    return _SS_error_flag;

}

_SS_ErrorFlag
PetscUI::InitVector( const KSP &A, std::unique_ptr<Vec> &x ) 
{

    PetscInt n,m;
    Mat mat_op, mat_pre;
    KSPGetOperators(A,&mat_op,&mat_pre);
    MatGetSize(mat_op, &m, &n);

    x.reset( new Vec() );
    VecCreate(MPI_COMM_WORLD, x.get() );
    VecSetSizes(*x, PETSC_DECIDE, m );
    VecSetFromOptions(*x);

    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::FreeVector( std::unique_ptr<Vec> &x) 
{
    VecDestroy(x.get());
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::FreeMatrix( std::unique_ptr<KSP> &A ) 
{
    KSPDestroy(A.get());
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::SetVector( std::unique_ptr<Vec> &x, std::vector<int> &cols, const std::string &type ) 
{
    /* check the vector size and cols are the same */
    PetscInt size;
    VecGetLocalSize( *x, &size );
    std::vector< double > vals(size,0);

    for ( int i = 0; i < size; i++ )
    {
        if ( cols [i] == 1 )
        {
            if ( type == "random" ) vals[i] = rand()/RAND_MAX;
            else if (type == "ones" ) vals[i] = 1.0;
        }
    }
    std::vector<int> cnums(size);
    std::iota(std::begin(cnums),std::end(cnums),0);
    VecSetValues(*x, size, cnums.data(),vals.data(), INSERT_VALUES);
    VecAssemblyBegin(*x);
    VecAssemblyEnd(*x);
    return _SS_error_flag;
}

_SS_ErrorFlag
PetscUI::Finalize( )
{
    KSPDestroy(&_ksp);
    return _SS_error_flag;
}

PetscErrorCode petscUI::KSPSetUp_SS(KSP ksp)
{

  PetscFunctionBegin; 
  PetscFunctionReturn(0);
}

/* This function is called when petsc calls "solve" */
PetscErrorCode petscUI::KSPSolve_SS(KSP ksp)
{
  PetscSolverSelecter *ss = ( PetscSolverSelecter* )ksp->data;     

  PetscFunctionBegin;
  Vec x,b;
  KSPGetSolution(ksp, &x );
  KSPGetRhs(ksp, &b );
  ss->Solve(ksp, x, b );

  PetscFunctionReturn(0);
}

/* Called when Petsc needs to destroy the ksp */
PetscErrorCode petscUI::KSPDestroy_SS(KSP ksp)
{
  PetscFunctionBegin;
  PetscSolverSelecter *ss = ( PetscSolverSelecter* )ksp->data;     

  ss->Finalize();
  delete ss; 
  PetscFunctionReturn(0);
}

/* Called when petsc wants to view the ksp solver */
PetscErrorCode petscUI::KSPView_SS(KSP ksp, PetscViewer viewer)
{
  PetscFunctionBegin;
  PetscSolverSelecter *ss = ( PetscSolverSelecter* )ksp->data;     

  ss->PrintChosenSolver();
  PetscFunctionReturn(0);
}

/* called when KSPSetFromOPtions is called inside Petsc */
PetscErrorCode petscUI::KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

/* Creates the petsc KSP solver with the solver selecter internals */
PetscErrorCode petscUI::KSPCreate_SS(KSP ksp)
{
  PetscErrorCode ierr;
  PetscSolverSelecter *ss;
  PetscFunctionBegin;

  ss = new PetscSolverSelecter();
  ksp->data = (void*)ss;

  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,3);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NATURAL,PC_LEFT,2);CHKERRQ(ierr);

  ksp->ops->setup          = petscUI::KSPSetUp_SS;
  ksp->ops->solve          = petscUI::KSPSolve_SS;
  ksp->ops->destroy        = petscUI::KSPDestroy_SS;
  ksp->ops->view           = petscUI::KSPView_SS;
  ksp->ops->setfromoptions = petscUI::KSPSetFromOptions_SS;
  ksp->ops->buildsolution  = KSPBuildSolutionDefault;
  ksp->ops->buildresidual  = KSPBuildResidualDefault;
  PetscFunctionReturn(0);
}


/* This is a static function that takes care of initializing the solver selecter as a petsc solver selecter */
_SS_ErrorFlag petscUI::InitializePetscSolverSelecter( KSP &ksp, std::shared_ptr<PetscUI> _interface )
{    
      /* Register the Petsc Solver Selecter */
      KSPRegister("KSPSS", petscUI::KSPCreate_SS);
      /* Set the type of the KSP to our SS type */
      KSPSetType(ksp, "KSPSS");      
      /* Get the preconditioner */
      PC pcc;
      KSPGetPC(ksp,&pcc);
      /* Turn of preconditioning in the solver selecter */
      PCSetType(pcc,PCNONE);

      /* Set up the KSP solver */
      KSPSetFromOptions(ksp); /* Not impliented yet */
      KSPSetSkipPCSetFromOptions(ksp,PETSC_TRUE);
      KSPSetUp(ksp); /* This builds the data structures */

      /* Get the data structure and use it to set the name of the database, then
       * to create and initialize the database. */
      PetscSolverSelecter *ss = ( PetscSolverSelecter* )ksp->data;     
      ss->Initialize(_interface);

      return 0;
}

#endif
