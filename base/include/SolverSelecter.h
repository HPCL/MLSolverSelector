#ifndef SS_SOLVERSELECTER_H
#define SS_SOLVERSELECTER_H

#include "ParameterBase.h"
#include "Solvers.h"
#include "Measurements.h"
#include "Features.h"
#include "GColoring.h"
#include "MachineLearning.h"
#include "DataBase.h"
#include "UserInterface.h"

/** \file SolverSelecter.h
 * \brief Header file for the solver selecter class
 */

/**
 * The SolverSelecter class is the main user implimented API for the entire project. The
 * user must impliment a derived class of the solver selecter with Matrix and vector set
 * to their problem specific matricies and vector data structures. For example, in petsc,
 * the call  is derived-class : public _SS_SolverSelecter<Mat,Vec>(...). The user must also
 * impliment some pure virtual functions, dicussed below
 * \tparam Matrix The matrix data structure to be used.
 * \tparam Vector The vector data structure to be used.
 **/
template <typename Matrix, typename Vector>
class _SS_SolverSelecter : public _SS_Parameters
{

private:

#if WITH_ZOLTAN
    //Only build a zoltan graph if compiled with Zoltan
    _SS_ZoltanGraphColor    *graph;    /**< Used to color a graph based on the sparcity (alpha) */
#endif

    std::shared_ptr<_SS_MachineLearning> machinemodel;                         /**< pointer the the machine learning class */
    std::shared_ptr< _SS_DataBaseBase >  database;                             /**< pointer to the database */
    _SS_Solver _solver;                                                        /**< store the solver for printing purposes */
    _SS_Measurements _measurements;                                            /**< measurements used by the solver selecter  */
    std::map<std::string,double> _mstruct;                                     /**< data required by the measurements  */
    _SS_Features _features;                                                    /**< features to be collected from each matrix  */
    std::string _mfile;                                                        /**< the current matrix filename   */
    std::set< std::string > filenames;                                         /**< list of matricies to test  */
    std::vector< std::string > solvers;                                        /**< list of solvers to test  */
    std::vector< std::string > preconds;                                       /**< list of preconditioners to test  */
    std::vector< std::map< std::string, std::set< std::string >> > parameters; /**< list of parameter sets to test   */
    std::unique_ptr<Matrix> A;                                                 /**< the current matrix */
    std::unique_ptr<Vector> x, b;                                              /**< the current vector */
    std::vector<int> ones;                                                     /**< a vector of ones for some reason */
    int nrows, ncols;                                                          /**< the number of rows and columns in the current matrix */
    std::vector< std::pair< unsigned int , unsigned int > > sparcity ;         /**< the sparcity pattern of the current matrix */
    std::vector< double > values ;                                             /**< the values associated with the current sparcity pattern */
    std::vector< std::pair < std::string, std::string > > sset ;               /**< set containing the allowed solver pairs */


public:
    std::shared_ptr<_SS_UserInterface<Matrix,Vector>> interface; /**< the user interface implimenting all this stuff */

    /** Constructor */
    _SS_SolverSelecter();

    /** Destructor */
    virtual ~_SS_SolverSelecter() ;

    /**  Initialize the solver selecter.  */
    virtual _SS_ErrorFlag Initialize(std::shared_ptr<_SS_UserInterface<Matrix,Vector>> _interface);

    /** Finalize the solver selecter  */
    virtual _SS_ErrorFlag Finalize();

    /** Turn solver selection on or off. If solver selection is turned off, the user is responsible
    *  for determining the best solvers to use. For example, the _SS_Solver returned in ExtractFeatures
    *  will be empty. Set using SetSolverSelection(bool) function. Note that ML should not be turned of
    *  when using the database builder.
    **/
    _SS_ErrorFlag SetSolverSelection(bool on);
    /** This is a conveinence based feature. Since the solver selecter has direct acess to the assembled
     *  matrix, it is very easy to equip with the ability to dump matricies. To use this feature, the user
     *  must impliment the virtual DumpMatrix() function. This parameter is set using the "SetMatrixDump()
     *  command
     **/
    _SS_ErrorFlag SetMatrixDump(bool on);


    /** Get features for a matrix free matrix */
    _SS_ErrorFlag ExtractFeatures( Matrix &A, std::string &mfile );

    /**
     * This is a convience feature that allows the user to add rows to the database inline. If this is
     * set, the parameter space described in the input_file will be used to solve the system. To level the
     * playing field, we use a random initial guess, and a RHS of ones, to complete the test, before continuing
     * on to solve the correct problem. As such, to use this feature, the user must impliment several optional
     * interface functions regarding initialization of the vectors, and setting the values.
     **/
    _SS_ErrorFlag SetBuildDatabaseInline( bool on,  std::string input_file );

    /** This is a convienece feature that allows the user to add rows to the database using an input file
     * and a list of saved matricies. In this case, the user must supply a function that opens the matrix
     * from file and initializes it. Given this, we can recurse over the parameter space listed in the input
     * file, makeing measurements as we go. ( This is what used to be the "DataBaseBuilder" );
     **/
    _SS_ErrorFlag BuildDataBaseFromFile( const std::string &inputfile );

    /** Return a pointer to the database controlled by the solver selector */
    _SS_ErrorFlag GetDataBase( _SS_DataBaseBase **_dat );


    /** Return a pointer to the interface */
    std::shared_ptr<_SS_UserInterface<Matrix,Vector>> GetInterface();



    /** Get the machine learner. Useful for setting ML parameters */
    _SS_ErrorFlag GetMachineLearner( _SS_MachineLearning *machine ) const;

    /** Add a solver to the list of solvers allowed to be selected by the solver selecter */
    _SS_ErrorFlag AddAllowedSolver( const std::string &solver, const std::string &precond );

    /** Print the parameter set for the chosen solver */
    _SS_ErrorFlag PrintChosenSolver();

    /** Select solver and solve the linear system. This is the main public function used
     * to solve the system. Inside here, we do the classification, ml, etc  */
    _SS_ErrorFlag Solve( Matrix &AA, Vector &xx, Vector &bb );


private:

    /** Train a model based on the data base */
    _SS_ErrorFlag Train();

    /**
    * recurse over the parameter space.
    **/
    _SS_ErrorFlag RecurseParameterSpace(Matrix &AA,  std::map< std::string, std::set< std::string >> &solver_params );
    /**
     * write the output to the database
     **/
    _SS_ErrorFlag WriteToDataBase( );

    /**
     * Classify a solver as good as bad
     **/
    _SS_ErrorFlag ClassifySolvers( );
    /**
     * Parse the parameter space
     **/

    _SS_ErrorFlag ParameterSpaceFromFile(std::string fname );


    /**
     * Add a matrix to list of matricies to test
     **/
    _SS_ErrorFlag AddMatrix( std::string filename );


    /**
     * Add a solver to the list of tests
     **/
    _SS_ErrorFlag AddSolver( std::string solver,
                             std::string preconditioner,
                             std::map< std::string, std::set< std::string >> &sparameters,
                             std::map< std::string, std::set< std::string >> &pparameters );

    /**
     * Split the strings based on a delimiter
     **/
    _SS_ErrorFlag StringSplit( const std::string &s, const char *delim, std::vector< std::string > &result );

};

#include "SolverSelecter.T"


#endif
