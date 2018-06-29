#ifndef DEFAULTPETSCFEATURE_SET
#define DEFAULTPETSCFEATURE_SET

#ifdef WITH_PETSCUI 

#include "typedefs.h"
#include "petscksp.h"
#include "petsc/private/kspimpl.h"
#include "petscdmda.h"

#define FULL

// Define a tolerance for a value being non-zero
#define NONZEROTOLERANCE 1e-125

// Define a processor to do the final calcuations on
#define ROOTPROC 0
#define HAVE_DIAGONAL 0

// Define the points for the samples.
#define POINTS_LOW 2
#define POINTS_HIGH 3
#define POINTS_LOW1 5
#define POINTS_HIGH1 7
#define POINTS_OFFSET 5


#ifdef FULL

#define DIMENSION 1
#define NNZ 1
#define AVGNONZEROSPERROW 1
#define ABSOLUTENONZEROSUM 1
#define FROBENIUSNORM 1
#define AVERAGEDIAGONALDISTANCE 2
#define ROWVARIANCE 1
#define TRACE 1
#define ABSOLUTETRACE 1
#define DIAGONALMEAN 1
#define DIAGONALAVERAGE ( 1 )
#define DIAGONALNONZEROS ( 1 || DIAGONALAVERAGE )
#define INFINITYNORM 1
#define MAXNONZEROSPERROW 1
#define DIAGONALSIGN 3
#define LOWERBANDWIDTH 1
#define UPPERBANDWIDTH 1
#define ROWDIAGONALDOMINANCE 1
#define MINNONZEROSPERROW 1
#define ONENORM 1
#define COLDIAGONALDOMINANCE 1
#define COLUMNVARIANCE 2
#define SYMMETRICITY 1
#define NONZEROPATTERNSYMMETRY 1
#define SYMMETRICINFINITYNORM 1
#define SYMMETRICFROBENIUSNORM 1
#define ANTISYMMETRICINFINITYNORM 1
#define ANTISYMMETRICFROBENIUSNORM 1

#endif

#ifdef RS1
#define DIMENSION 1
#define NNZ 1
#define AVGNONZEROSPERROW 0
#define ABSOLUTENONZEROSUM 0
#define FROBENIUSNORM 0
#define AVERAGEDIAGONALDISTANCE 0
#define ROWVARIANCE 0
#define TRACE 0
#define ABSOLUTETRACE 1
#define DIAGONALMEAN 0
#define DIAGONALAVERAGE ( 1 )
#define DIAGONALNONZEROS ( 1 || DIAGONALAVERAGE )
#define INFINITYNORM 1
#define MAXNONZEROSPERROW 0
#define DIAGONALSIGN 0
#define LOWERBANDWIDTH 1
#define UPPERBANDWIDTH 0
#define ROWDIAGONALDOMINANCE 0
#define MINNONZEROSPERROW 1
#define ONENORM 0
#define COLDIAGONALDOMINANCE 0
#define COLUMNVARIANCE 2
#define SYMMETRICITY 0
#define NONZEROPATTERNSYMMETRY 1
#define SYMMETRICINFINITYNORM 0
#define SYMMETRICFROBENIUSNORM 0
#define ANTISYMMETRICINFINITYNORM 0
#define ANTISYMMETRICFROBENIUSNORM 0
#endif

#ifdef RS2
#define DIMENSION 1
#define NNZ 1
#define AVGNONZEROSPERROW 0
#define ABSOLUTENONZEROSUM 0
#define FROBENIUSNORM 0
#define AVERAGEDIAGONALDISTANCE 0
#define ROWVARIANCE 0
#define TRACE 0
#define ABSOLUTETRACE 1
#define DIAGONALMEAN 0
#define DIAGONALAVERAGE ( 1 )
#define DIAGONALNONZEROS ( 1 || DIAGONALAVERAGE )
#define INFINITYNORM 0
#define MAXNONZEROSPERROW 0
#define DIAGONALSIGN 0
#define LOWERBANDWIDTH 1
#define UPPERBANDWIDTH 0
#define ROWDIAGONALDOMINANCE 0
#define MINNONZEROSPERROW 1
#define ONENORM 0
#define COLDIAGONALDOMINANCE 0
#define COLUMNVARIANCE 2
#define SYMMETRICITY 0
#define NONZEROPATTERNSYMMETRY 0
#define SYMMETRICINFINITYNORM 0
#define SYMMETRICFROBENIUSNORM 0
#define ANTISYMMETRICINFINITYNORM 0
#define ANTISYMMETRICFROBENIUSNORM 0
#endif

#define SUMCOUNT (NNZ + AVGNONZEROSPERROW  + ABSOLUTENONZEROSUM + FROBENIUSNORM  + \
  AVERAGEDIAGONALDISTANCE + ROWVARIANCE + TRACE + ABSOLUTETRACE + DIAGONALNONZEROS +DIAGONALMEAN +  DIAGONALAVERAGE )


#define DIAGONAL (TRACE || ABSOLUTETRACE || DIAGONALMEAN || DIAGONALSIGN || DIAGONALNONZEROS || DIAGONALAVERAGE)

#define MAXCOUNT (INFINITYNORM + MAXNONZEROSPERROW + DIAGONALSIGN + LOWERBANDWIDTH + UPPERBANDWIDTH + ROWDIAGONALDOMINANCE )

#define MINCOUNT (MINNONZEROSPERROW)

#define COLCOUNT (ONENORM + COLDIAGONALDOMINANCE + COLUMNVARIANCE )

#define SYMMETRY (SYMMETRICITY || NONZEROPATTERNSYMMETRY || SYMMETRICFROBENIUSNORM || \
        SYMMETRICINFINITYNORM || ANTISYMMETRICINFINITYNORM || ANTISYMMETRICFROBENIUSNORM)

#define FEATURECOUNT(npoints) MAXCOUNT + MINCOUNT + COLCOUNT*npoints + SUMCOUNT + SYMMETRY*npoints*npoints

namespace SolverSelecter
{


class PetscTestingSpace
{
public:

    PetscTestingSpace( );

    virtual void
    extract_features(KSP &ksp,
                     std::map<std::string, double> &fmap ) = 0;

    virtual void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) = 0;

    virtual void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double > &mmap) = 0;

    virtual void
    classify_measurements(std::map<std::string, double> &cmap ) = 0;

};


class DefaultPetscTestingSpace : public PetscTestingSpace
{
public:

    std::chrono::time_point<std::chrono::high_resolution_clock> time_start;

    void
    extract_features( KSP &ksp,
                      std::map<std::string, double> &fmap ) override;

    void
    start_measurements(KSP &ksp,
                       Vec &x,
                       Vec &b ) override;

    void
    stop_measurements(KSP &ksp,
                      Vec &x,
                      Vec &b,
                      std::map<std::string, double> &mmap ) override;

    void
    classify_measurements(std::map<std::string, double> &cmap ) override;

    int
    MapDMNumberingToNaturalNumbering(std::vector< std::pair< int , int > > &points );

    int
    GetSamplePoints(int n,
                    int edge,
                    int interior,
                    std::vector<std::pair<int,int>> &points );

    int
    GetJacobianColumn(Mat J,
                      std::pair< int, int > &point,
                      Vec *s );

    int
    GetJacobianSamplePoints(Mat J,
                            int edge,
                            int interior,
                            std::vector< std::pair< int, int > >&points );

    static void
    MPI_FeatureReduce(void *invec,
                      void *inoutvec,
                      int *len,
                      MPI_Datatype *datatype);

    int
    ExtractJacobianFeatures(Mat J,
                            int edge,
                            int interior,
                            std::map< std::string, double > &fnames );

};

}
#endif

#endif
