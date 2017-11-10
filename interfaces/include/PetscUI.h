#ifndef PETSCUI_HEADER
#define PETSCUI_HEADER

/** \file PetscSolverSelecter.h
 * \brief In this example, we precompile a solver selecter based on Petsc. This creates a library
 * that can be used in other example. This is different from the dummy example, where the entire 
 * Solver selecter library is built at compilation time of the example.  
 **/ 

#ifdef WITH_PETSCUI

//Include the petsc files
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "SolverSelecter.h"
#include "MatrixNorm.h"
#include "MatrixVariance.h"
#include "CPUTime.h"
#include "Sqlite3.h"
#include "MLWaff.h"

namespace petscUI {


/* The first step is to define the PetscSolverSelecterClass. This is used in
 * the PetscUI defined below, so it needs to come first. For petsc, we use the 
 * ksp as the matrix and a petsc Vec as the vector. This gives us more data in the solve
 * which is nice. */

typedef _SS_SolverSelecter<KSP,Vec> PetscSolverSelecter;

/*Next, we need to impliment a UserInterface for the solver selecter. This is the main
 * component required by the solver selecter package.  
 */

class PetscUI : public _SS_UserInterface<KSP,Vec>
{
public :

  KSP _ksp; /**< linear solver */
  std::string inputfile;  /**< input file to set the solvers and preconditioners and matrix files */
  std::string dname; /**< the name of the database to use */ 
  bool solver_on; /**< should we use the solver selecter ? */

  bool dump_on; /**< should we dump the matricies to file */
  bool inline_on; /**< should we perform database building inline (i.e, runtime database build )*/

  PetscUI(std::string _db /**< database name */, 
          std::string _in,/**< input file */ 
          bool _so, /**< use the solver ?*/
          bool _do, /**< dump matricies ? */
          bool _io /**< build database during runtime based on input file */);
  
  /* The interface functions that must be overriden. Some are optional, but, we impliment them all */ 
  _SS_ErrorFlag AddFeaturesAndMeasurements( _SS_Features &features, _SS_Measurements &measure ) override;
  _SS_ErrorFlag GetDataBaseImpl( std::shared_ptr< _SS_DataBaseBase > &database );
  _SS_ErrorFlag GetMLImpl( std::shared_ptr< _SS_MachineLearning > &machinemodel );
  _SS_ErrorFlag SetSolverSelecterOptions( _SS_SolverSelecter<KSP,Vec> *ss ) override;  
  _SS_ErrorFlag GetMatrixInfo( KSP &A, int &nrows, int &ncols, int &chunks, std::string &matrix_name, bool &mfree ) override;    
  _SS_ErrorFlag MatVecAndGetData( KSP &ksp, Vec &x, std::vector<int> &row, std::vector<double> &val) override ;
  _SS_ErrorFlag SolveSystem( KSP &ksp, Vec &x, Vec &b, _SS_Solver &solver, std::map<std::string,double> &mstruct) override;  
  _SS_ErrorFlag InitMatrix( std::string filename, std::unique_ptr<KSP> &A ) override;
  _SS_ErrorFlag InitVector( const KSP &A, std::unique_ptr<Vec> &x ) override; 
  _SS_ErrorFlag FreeVector( std::unique_ptr<Vec> &x) override;
  _SS_ErrorFlag FreeMatrix( std::unique_ptr<KSP> &A ) override;
  _SS_ErrorFlag SetVector( std::unique_ptr<Vec> &x, std::vector<int> &cols, const std::string &type ) override; 
  _SS_ErrorFlag GetSparcity( KSP &ksp, int &chunk, std::vector< std::pair< unsigned int, unsigned int > > &sparcity,
                            std::vector< double > &values ) override ;


  _SS_ErrorFlag Finalize(); 

};



PetscErrorCode KSPSetUp_SS(KSP ksp);
PetscErrorCode KSPSolve_SS(KSP ksp);
PetscErrorCode KSPDestroy_SS(KSP ksp);
PetscErrorCode KSPView_SS(KSP ksp, PetscViewer viewer);
PetscErrorCode KSPSetFromOptions_SS(PetscOptionItems *PetscOptionsObject,KSP ksp);
PetscErrorCode KSPCreate_SS(KSP ksp);

_SS_ErrorFlag InitializePetscSolverSelecter( KSP &ksp, std::shared_ptr<PetscUI> _interface );

}

#endif
#endif
