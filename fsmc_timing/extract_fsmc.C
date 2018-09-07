
//


//


#include <sys/stat.h>
#include <unistd.h>
#include <fstream>      // std::fstream
#include <chrono>

#include "FeatureExtract.h"

static char help[] = "Loads a Petsc Matrix in PetscBinary form and extracts the features.\n\n";

int main(int argc,char **argv)
{

  PetscInitialize(&argc,&argv,(char*)0,help);

 
  if ( argc < 5 ) {
    std::cout << " Useage " << argv[0] << " <path-to-matrix> <output_file> <edgepoints> <interiorpoints> " << std::endl;
    return 1;
  }

  std::string matrix_file = argv[1];
  std::string output_file = argv[2];
  
  int edge = std::atoi(argv[3]);
  
  double interior = std::atof(argv[4]);
  bool percentage = false;

  auto start = std::chrono::high_resolution_clock::now();
 
  // Load the Matrix 
  Mat A;
  PetscViewer fd;
  MatCreate(PETSC_COMM_WORLD, &A);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_file.c_str(), FILE_MODE_READ, &fd);
  MatSetFromOptions(A);
  MatLoad(A,fd);
  PetscViewerDestroy(&fd);

  Vec x,b;
  MatCreateVecs(A, &x, &b );
  VecSet(b,1.0) ;
  VecSetRandom(x,NULL);
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  
  auto stop1 = std::chrono::high_resolution_clock::now();
  
  KSPSetType(ksp, "gmres");
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, "jacobi");
  KSPSetOperators(ksp,A,A);
  
  auto stop2 = std::chrono::high_resolution_clock::now();

  //Extract the features
  std::vector< std::pair<std::string, double> > feature_set;
  
  
  int success = ExtractJacobianFeatures( A, edge, interior, feature_set, 1 );


  auto stop3 = std::chrono::high_resolution_clock::now();
 
  long long mload    = std::chrono::duration_cast<std::chrono::microseconds>(stop1-start).count(); 
  long long kspsetup = std::chrono::duration_cast<std::chrono::microseconds>(stop2-stop1).count(); 
  long long extract  = std::chrono::duration_cast<std::chrono::microseconds>(stop3-stop2).count(); 
   
  //Destroy the matrix
  MatDestroy(&A);

  // THis is just data output junk   
  int size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank == 0 ) {

    std::ostringstream feature_labels;
    feature_labels << "Matrix Name, Comm Size , Edge, Interior "; 

    std::ostringstream feature_values;
    feature_values << matrix_file;
    feature_values << "," << size;
    feature_values << "," << edge;
    feature_values << "," << interior;
     
      
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
   

    std::fstream fs, fss;
    fs.open(output_file, std::fstream::out | std::fstream::app ) ;
    fss.open(output_file + "_timing", std::fstream::out | std::fstream::app ) ;
    
    if ( ! exists ) {
      fs << feature_labels.str(); 
      fss << " Matrix , Load , Setup , Extract  \n " ;
    }

    fs << feature_values.str() ; 
    
  std::ostringstream tim; 
    tim << matrix_file << "," << mload << "," << kspsetup << "," << extract << "\n";
    fss << tim.str();
    fss.close(); 
    fs.close();
  }

  PetscFinalize() ;
  return 0;
}
