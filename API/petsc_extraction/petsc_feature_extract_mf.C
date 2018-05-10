#include "PetscFeatureExtraction.h"

//Currently, we get different results for the cases where the matrix
//is blocked vs striped. This 
#define BLOCKED PETSC_DECIDE // 1 toggle u


/** \file ex1.C
 * \brief Extract features from a matrix free petsc shell instance. where F(u(x,y)) = Soutce **/

static char help[] = "Loads a Petsc Matrix in PetscBinary form and extracts the features.\n\n";


typedef struct {
  PetscInt mx = 5 ;   
  PetscInt my = 5 ; 
  void *appctx; 
  
  // current is the vector about which we are estimating the derivitive
  // current calc is A(current). This is used repeatedly, so we save it. 
  Vec current, current_calc;  

  // Function pointers to set up the problem. 
  void (*Evaluate)( Mat mat, Vec x, Vec y ) = NULL; 
  double (*Source)(double u, void *ctx ) = NULL; // Provide the source term at the current postion 
  double (*Boundary)( void ) = NULL; // Only support constant Dirch boundary terms for now. 

} MFMatrix;

int destroy_mf_matrix( MFMatrix * appctx )
{
   VecDestroy( &(appctx->current) );
   VecDestroy( &(appctx->current_calc) );
}


// Set the Vector about which we will take the derivitive
int SetCurrentState( Mat mat, Vec x ) 
{
  void * ctx;
  MatShellGetContext( mat, &ctx );
  MFMatrix *appctx = ( MFMatrix * ) ctx ;
  
  VecDuplicate(x, &(appctx->current));
  VecDuplicate(x, &(appctx->current_calc));
  VecCopy(x,appctx->current); 

  //Evaulate it now because it will be used a lot later. 
  appctx->Evaluate(mat,appctx->current,appctx->current_calc); 
}


//This is the function that gets called whenever a MatVec is needed. 
int MatrixFreeJacobian( Mat A, Vec x, Vec y  ) // Ax = y 
{
    void * ctx;
    MatShellGetContext( A, &ctx );
    MFMatrix *appctx = ( MFMatrix * ) ctx ;
    
    //Calcuate the perturbation factor
    PetscInt n ;
    VecGetSize( x , &n ) ; 
    PetscScalar per =  sqrt(1e-16)*2.0/ (double) n ;
    
    //Calcuate the preturbed vector (also learn to spell)
    Vec preterbed; 
    VecDuplicate( x, &preterbed );
    VecWAXPY( preterbed, per, x, appctx->current );      // preturbed = ( current + per * x ) 
    
    //Evaulate preturbed vector
    appctx->Evaluate( A, preterbed, y );  // y = F( current + per * x ) 
    
    //Calcuate the jacobian difference 
    VecAXPY( y , -1.0, appctx->current_calc ) ;   //Ax = y = F( current + a x ) - F(current)  
    VecScale( y , 1.0/per );
    VecDestroy(&preterbed);
    return 0;
}

int create_shell_matrix( MFMatrix *appctx, DM *dm, Mat *Jaco, Vec *about )
{

   PetscInt m = appctx->mx;
   PetscInt n = appctx->my;
   
   //Create a 2d petsc DM. BLOCKED is a Macro that sets the number of processors in the x direction. 
   DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,appctx->mx,appctx->my,
                BLOCKED,PETSC_DECIDE,1,1,0,0,dm);
   
   //Set up the DM for a unit box.
   DMSetFromOptions(*dm);
   DMSetUp(*dm);
   DMDASetUniformCoordinates(*dm,0.0,1.0,0.0,1.0,0.0,1.0);   
   
   //Build the shell matrix and set the matrix-vec function . 
   DMSetMatType(*dm, MATSHELL );
   DMCreateMatrix(*dm, Jaco);
   MatShellSetContext(*Jaco, appctx );
   MatShellSetOperation(*Jaco, MATOP_MULT, (void (*)(void))MatrixFreeJacobian);
   
   //Create a global vector for this DM. 
   DMCreateGlobalVector(*dm, about);
}

int destroy_shell_matrix( MFMatrix *appctx, DM *dm, Mat *Jaco, Vec *about )
{
  //Destroy shell matrix and ssociated structures. 
  VecDestroy(about);
  MatDestroy(Jaco);
  DMDestroy(dm);
}

 
// Impls for the bratu problem source term  
double Bratu_Source(double u, void *ctx )
{       
    double *lambda = (double*) ctx;
    return exp(u)*(*lambda); 
}

// The boundary is zero. 
double DBoundary(void)
{
  return 0.0;
}

//Calculate y = F(x) where F  = uxx + uyy + source = 0 
void Poisson(Mat mat, Vec x, Vec y )
{
    void * ctx;
    MatShellGetContext( mat, &ctx );
    MFMatrix *appctx = ( MFMatrix * ) ctx ;
    
    DM dm;
    MatGetDM( mat, &dm );
   
    PetscInt Mx, My;
    DMDAGetInfo(dm,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                   PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);

    PetscScalar hx     = 1.0/(PetscReal)(Mx-1);
    PetscScalar hy     = 1.0/(PetscReal)(My-1);
    PetscScalar hxdhy  = hx/hy;
    PetscScalar hydhx  = hy/hx;
    
    // Create a local vector and populate its ghost cells. 
    Vec local_x, local_y; 
    DMCreateLocalVector(dm, &local_x );
    DMGlobalToLocalBegin(dm, x, INSERT_VALUES, local_x);
    DMGlobalToLocalEnd(dm, x, INSERT_VALUES, local_x);
   
    // Get a local version to store the solution -- set to zero 
    DMCreateLocalVector(dm, &local_y );
    VecSet(local_y,0.0);
    
    // Get the underlying data arrays for the local vectors. DMDAVecGetArray is 
    // nice because it extracts a 2D array.  
    PetscScalar **xarray, **yarray ;
    DMDAVecGetArray(dm, local_x, &xarray);
    DMDAVecGetArray( dm, local_y, &yarray ); 

    PetscInt xs,ys,xm,ym, rank ;
    DMDAGetCorners(dm,&xs,&ys,NULL,&xm,&ym,NULL);
    PetscScalar u,un,us,ue,uw,uxx,uyy;
    for ( int j = ys; j < ys + ym; j++ ) {
      for ( int i = xs; i < xs + xm; i++ ) {
        if ( i == 0 || j == 0 || i == Mx-1 || j == My-1 ) yarray[j][i] = appctx->Boundary(); 
        else {
          u = xarray[j][i] ; 
          un = xarray[j-1][i];
          us = xarray[j+1][i];
          ue = xarray[j][i-1];
          uw = xarray[j][i+1];
          uxx = (2.0*u - ue - uw )*hydhx ;
          uyy = (2.0*u - un - us )*hxdhy ;
          yarray[j][i] = uxx + uyy ;
          if ( appctx->Source != NULL ) 
              yarray[j][i] -= hx*hy*appctx->Source(u, appctx->appctx ); 
              
        }
      }
    }
    
    DMDAVecRestoreArray(dm, local_x, &xarray );
    DMDAVecRestoreArray(dm, local_y, &yarray );
    
    //Distribute the local vector to the global vectors
    VecSet(y, 0.0);
    DMLocalToGlobalBegin(dm, local_y, ADD_VALUES, y );
    DMLocalToGlobalEnd(dm, local_y, ADD_VALUES, y );

    //Clean up
    VecDestroy(&local_x);
    VecDestroy(&local_y);
}


/** Main */
int main(int argc,char **args)
{

  PetscInitialize(&argc,&args,(char*)0,help);
  
  //do we want to add the bratu source term  
  bool bratu = true;
  
  //Create our matrix context 
  MFMatrix appctx; 
  appctx.mx = 5;
  appctx.my = 5;
  appctx.Evaluate = &Poisson; 
  appctx.Boundary = &DBoundary;
  
  // Put the bratu parameter into the context. 
  double lambda = 2.5;
  if ( bratu  )
  {    
      appctx.Source = &Bratu_Source;
      appctx.appctx = (void*) &lambda; 
  }

  //Create the DM, Mat and Vec 
  DM dm;
  Mat Jaco;
  Vec about;
  create_shell_matrix( &appctx , &dm, &Jaco, &about );
  
  //Choose a vector to complete the derivitive about. In this case a vector of ones. Note that
  //this vector doesn't conform to the boundary condition. The Poisson function deals with this,
  //but the matlab version of the bratu problem does not. 
  VecSet(about, 1.0 ); 
  SetCurrentState(Jaco, about );
  
  //Extract the features.  
  ExtractJacobianFeatures(Jaco);
  
  destroy_shell_matrix( &appctx, &dm, &Jaco, &about );
  destroy_mf_matrix(&appctx);

  PetscFinalize() ;
  return 0;
}



