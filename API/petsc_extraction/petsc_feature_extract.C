#include "PetscFeatureExtraction.h"


static char help[] = "Loads a Petsc Matrix in PetscBinary form and extracts the features.\n\n";

int main(int argc,char **args)
{

  PetscInitialize(&argc,&args,(char*)0,help);
 
  // Matrix file name   
  std::string filename("./data/662_bus.pbin");
 
  // Load the Matrix 
  Mat A;
  PetscViewer fd;
  MatCreate(PETSC_COMM_WORLD, &A);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &fd);
  MatSetFromOptions(A);
  MatLoad(A,fd);
  PetscViewerDestroy(&fd);

  //Extract the features
  ExtractJacobianFeatures( A );
  
  //Destroy the matrix
  MatDestroy(&A);

  PetscFinalize() ;
  return 0;
}
