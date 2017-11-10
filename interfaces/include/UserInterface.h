#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Features.h"
#include "Measurements.h"
#include "Solvers.h"
#include "DataBase.h"
#include "MachineLearning.h"


template<typename Matrix, typename Vector>
class _SS_SolverSelecter;

template<typename Matrix,typename Vector>
class _SS_UserInterface
{
public:
    /* This is a virtual base class that must be overwritten
     * by the user of the code. It holds all the user defined functions
     * required by the solver selecter. */
    _SS_UserInterface();

    ///////////// Solver Selecter Options and settings /////////////////////////

    /** Add features that will be tested */
    virtual _SS_ErrorFlag AddFeaturesAndMeasurements( _SS_Features &features, _SS_Measurements &measure );

    virtual _SS_ErrorFlag SetSolverSelecterOptions( _SS_SolverSelecter<Matrix,Vector> *ss );

    /** User must choose a database implimentation. The sqlite3 implimentation is provided */
    virtual _SS_ErrorFlag GetDataBaseImpl( std::shared_ptr< _SS_DataBaseBase > &database ) = 0;

    /** The user must also set the ML implimentation */
    virtual _SS_ErrorFlag GetMLImpl( std::shared_ptr< _SS_MachineLearning > &machinemodel /**< ptr to init*/ ) = 0;

    /** User impliments this solve function. The options chosen by the solver
    * selecter can be obtained by calling solvers.GetChosenOptions(...)  */
    virtual _SS_ErrorFlag SolveSystem( Matrix &A, Vector &x, Vector &b,
                                       _SS_Solver &solver,
                                       std::map< std::string, double > &mstruct ) = 0 ;

    /** virtual function implimented by the user to get the matrix information. Chunck should be set to minimize
     * the memory consumption. Basically, the user must return a vector indicating the sparicity structure
     * for the matrix. If the matrix is huge, we cannot allocate enough space for the full matrix, so the sparcity
     * pattern is read in chuncks, before being discarded. TODO Slow??? */
    virtual _SS_ErrorFlag GetMatrixInfo( Matrix &A, int &nrows, int &ncols, int &chunks, std::string &matrix_name, bool &mfree );

    /** virtual function implimented by the user to get the sparcity pattern. Chunk indicates that chunk of data
     * that we are requesting as specified in the SolveSystem call. If the matrix is matric free, the sparcity
     * pattern is all that is required. */
    virtual _SS_ErrorFlag GetSparcity( Matrix &A, int &chunk,
                                       std::vector< std::pair< unsigned int, unsigned int > > &sparcity,
                                       std::vector< double > &value );

    /** Only needed if the Zoltan graph partitioner is used to fill in matrix values for matrix free
     * feature extraction. **/
    virtual _SS_ErrorFlag MatVecAndGetData( Matrix &A, Vector &x, std::vector<int> &rows, std::vector<double> &vals );

    /** (optional) dump matrix to file **/
    virtual _SS_ErrorFlag DumpMatrix( Matrix &A );

    /**
     *  This is an optional function only required if using the BuildDataBaseFromFile option.
     * \param op_filename Filename of the operator matrix to open up
     **/
    virtual _SS_ErrorFlag InitMatrix( std::string filename,
                                      std::unique_ptr<Matrix> &A );

    /** Initialize a vector based on the Matrix in AA */
    virtual _SS_ErrorFlag InitVector( const Matrix &AA,
                                      std::unique_ptr<Vector> &x );

    /**
     * Free the pointers needed during building of the database.
     **/
    virtual _SS_ErrorFlag FreeMatrix( std::unique_ptr<Matrix> &A );

    /** Override if some actions must be completed before the vector is destroyed. I.e, for Petsc,
     * we must vall VecDestroy before the pointer can be destroyed. */
    virtual _SS_ErrorFlag FreeVector( std::unique_ptr<Vector> &x );

    /**
     * (optional) Set the vector based on the matrix. cols is a vector listing the indices that need to be set,
     * and type is what they should be set to ( "random", "ones", or "zero" );
     TODO make cols a const vector, otherwise could get slammed by the user
     **/
    virtual _SS_ErrorFlag SetVector( std::unique_ptr<Vector> &cloned_x, std::vector<int> &cols, const std::string &type );

};

#include "UserInterface.T"

#endif

