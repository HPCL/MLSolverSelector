
/** \file ex01.C
 * \brief SS-LITE example. SSLITE is a stripped down version of the SolverSelector API. It allows you to 
 * use a pre-built C50 machine learning model to determine the best linear solver inside a petsc solve. 
 *
 * Using the library in a Petsc executable . 
 *
 * 1) Build the SS-Lite library ( see Readme in previous directory ) 
 * 2) Add a -I <ss-lite-install-dir/include  to your compile flags
 * 3) Add a -L <ss-lite-install-dir>/lib -lsslite to you library flags
 * 
 * 4) In your executaable add a '  #include "PetscInterface.h" ' at the top
 * 5) Somewhere after PetscInitialize is called : call the petsc coupler to register the ksp
 *
 *      PetscCoupler::CreateSolverSelectorKSP();
 *
 * 6) At this point, the special solver selector KSP has been linked with Petsc. 
 * 7) Using the ksp at run time is achieved by setting the appropriate command line, 
 *
 *    ./<executable>.exe ... -ksp_type KSPSS -pc_type pcnone -ksp_ss_filestem <C50 filestem> -ksp_ss_featureset <featureSetName> 
 *
 *    Some additional options are:
 *        -ksp_ss_edge <int>    --- This is the number of edge points to sample when extracting features
 *        -ksp_ss_interior <int> --- This is the number of interior points to sample when extracting features
 *        -ksp_ss_matvecs <true|false> --- Set this false to use memory access instead of matvecs to extract the sample columns. 
 *
 *   NOTES:
 *       a) the filestem and featureset parameters are required -- The program will abort if you choose the KSPSS solver but do not provide them
 *       b) The filestem represents the standard C50 filestem found from a prebuilt model. The files required are  
 *             i) <filestem>.data ( i don't actually know if we need this one )
 *             ii) <filestem>.names
 *             iii) <filestem>.tree 
 *             iv) IMPORTANT <filestem>.solvers
 *             
 *                    This is the file that links the hashes in the tree to the solver parameters required in petsc to set up 
 *                    the solver -- so it is essential that it is provided. The expected format is :
 *                           <hash> , <solver-name> <preconditioner-name> <parameter> <value> <parameter> <value> ...
 *                           <hash> , <solver-name> <preconditioner-name> <parameter> <value> <parameter> <value> ...
 *                            .
 *                            .
 *                            .
 *                     
 *                     Note: If you have a command line parameter that does not accept a value, then, you should
 *                     put -ksp_<parameter> SSNONE     
 *             
 *        c) The feature set value is also very important. This value tells sslite what feature set to use when extracting 
 *           values from the matrices. For the best accuracy with the lowest costs, it is important that the feature set 
 *           that is set using this parameter contains the same features used in the model. I.e., ithe column names in
 *           the <filestem>.names file must match those in the feature set (with no spelling mistakes). 
 *
 *           To set up a new feature set, see the README in the Feature sets directory. If there are some differences in the
 *           spelling of a feature in the model and the spelling of a feature in the code, or if a feature in the model 
 *           is not available in the code, an error should be thrown (todo). The use case will generally be that the 
 *           model is built using features extracted from the same feature set -- but in testing this will not always
 *           be the case. 
 *   
 *    The example below is a generic matrix solver found in the petsc examples. All modifications to the 
 *    file required to integrate SS-LITE are commented (there are two of them :) 
 *          
 *  ** FYI -- this example solves the matvec Ax = b where b is all ones for whatever
 *     matrix you pass in using the -f <matrix-file> parameter.  
 *
 *
 ******************** TODO ***************************************
 *  1) Put in code to ensure the matrix loaded from file has an explicit diagonal -- Some solvers in petsc need this. 
 *  2) Put in a MOOSE example. ( running moose out of source is strange?) 
 * 
 * **/
static char help[] = " Type -help to see the options. nano ex01.cpp to see a detailed readme.  \n \n";

/******** Here is the include required ********/
#include "PetscInterface.h"    // Normally there would be some other petsc includes here, but we already got them in this include. 
/*********************************************/



/** Main */
int main(int argc,char **args)
{

  PetscMPIInt    size;
  PetscErrorCode ierr; 
  PetscInitialize(&argc,&args,(char*)0,help);

  /****** Here is the call to the coupler. This regiesters the KSPSS solver for use  */
  PetscCoupler::CreateSolverSelectorKSP();
  /************ That is it, everything below here is plain old petsc.i Just make sure to set the 
   * appropriate command line parameters as described above.  ***********/ 
  
  char file[PETSC_MAX_PATH_LEN];
  PetscBool flg;
  PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,&flg);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f option");


  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  Vec            x, b, u;      /* approx solution, RHS, exact solution */
  PetscReal      norm;  /* norm of solution error */
  PetscScalar    neg_one      = -1.0,one = 1.0,value[3];
   
  Mat A;
  PetscViewer fd;
  MatCreate(PETSC_COMM_WORLD, &A);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file, FILE_MODE_READ, &fd);
  MatSetFromOptions(A);
  MatLoad(A,fd);
  PetscViewerDestroy(&fd);

  PetscInt n,m;
  MatGetSize(A, &m, &n);
  VecCreate(MPI_COMM_WORLD, &x );
  VecSetSizes(x, PETSC_DECIDE, m );
  VecSetFromOptions(x);  
  VecDuplicate(x,&b);
  VecDuplicate(x,&u);

  ierr = VecSet(u,one);CHKERRQ(ierr);
  ierr = MatMult(A,u,b);CHKERRQ(ierr);

  KSP ksp; 
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);  
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  
  //ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  //ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  
  PetscReal res, xnorm;
  MatMult(A,x,u);  // Ax = u;
  VecAXPY(u, neg_one, b ); // Ax -b = u 
  VecNorm(u,NORM_2,&res);
  VecNorm(x,NORM_2,&xnorm);
  double scaled_residual = (double)res / (double) xnorm ; 


  PetscInt its;  
  KSPGetIterationNumber(ksp,&its);
  PetscReal rnorm;
  KSPGetResidualNorm(ksp,&rnorm);

  PetscPrintf(PETSC_COMM_WORLD,"Final Relative Residual %g iterations %D\n",(double)scaled_residual,its);

  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return 0;
}
