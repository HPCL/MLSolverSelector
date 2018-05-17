#include "PetscFeatureExtraction.h"
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>      // std::fstream

static char help[] = "Loads a Petsc Matrix in PetscBinary form and extracts the features.\n\n";

int main(int argc,char **argv)
{

  PetscInitialize(&argc,&argv,(char*)0,help);

 
  if ( argc < 3 ) {
    std::cout << " Useage " << argv[0] << " <path-to-matrix> <output_file> " << std::endl;
    return 1;
  }

  std::string matrix_file = argv[1];
  std::string output_file = argv[2];
   
  // Load the Matrix 
  std::cout << " Loading the matrix " << matrix_file << std::endl;
  Mat A;
  PetscViewer fd;
  MatCreate(PETSC_COMM_WORLD, &A);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_file.c_str(), FILE_MODE_READ, &fd);
  MatSetFromOptions(A);
  MatLoad(A,fd);
  PetscViewerDestroy(&fd);

  std::cout << " Loaded the matrix " << matrix_file << std::endl;
  
  std::cout << " Starting feature extraction " << matrix_file << std::endl;
  //Extract the features
  std::vector< std::pair<std::string, double> > feature_set;
  ExtractJacobianFeatures( A, feature_set );
  

  std::cout << " Finished Feature extraction " << matrix_file << std::endl;

  //Destroy the matrix
  MatDestroy(&A);

  // THis is just data output junk   
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  std::ostringstream feature_labels;
  feature_labels << "Matrix, Comm Size "; 

  std::ostringstream feature_values;
  feature_values << matrix_file;
  feature_values << "," << size; 
    
  for ( auto it : feature_set ) {
      feature_labels << "," + it.first;
      feature_values << "," << it.second; 
  }
  feature_labels << "\n";
  feature_values << "\n";

  bool exists = false;
  if ( ( access( output_file.c_str(), F_OK ) != -1 ) ) 
  {
      exists = true; 
  }
 
  std::fstream fs;
  fs.open(output_file, std::fstream::out | std::fstream::app ) ;
  if ( ! exists ) {
    fs << feature_labels.str(); 
  }
  fs << feature_values.str() ; 
  fs.close();


  PetscFinalize() ;
  return 0;
}
