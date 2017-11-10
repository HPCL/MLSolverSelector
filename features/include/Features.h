#ifndef _SS_FEATURES_H
#define _SS_FEATURES_H

#include "ParameterBase.h"
/**
 * \file Features.h
 * \brief Header file for a "feature". The _SS_Feature class is a structure that allows
 * us to extract features from the matrix, independent of the implimentation.
 **/


/** Base class for a feature . A feature contains 3 key functions that
 * must be implimented; initialize, finalize and loop. FIXME This is pretty slow, might not
 * be the best way to go about it. Maybe better to leave feature imp in the interface */
class _SS_Feature
{
public:

    std::string name;
    /** Construnctor */
    _SS_Feature(const std::string &name) ;
    /* Initialize feature collection. This makes sure the feature is in the correct state
     * before calling the subclass initialization routine. */
    _SS_ErrorFlag Initialize( MPI_Comm comm /**< MPI communicator to be used with this feature */,
                              const int &num_rows /**< number of rows in the matrix */,
                              const int &num_cols /**< number of columns in the matrix */);

    /**< Next is the main function of the feature. Next gives a (row,column,data). This function
     * is called once for every entry in the sparcity pattern of the matrix. FIXME Might need to add a
     * flag for features that need to know all values rather than just non-zeros?
     **/
    _SS_ErrorFlag Next( MPI_Comm comm /**< communicator for this feature*/,
                        const int &row /**< row index */,
                        const int &col /**< col index */,
                        const double &data /**< data value */ );

    /** Finalize the feature collection. This is where MPI communication will be done to collect
     * the global value.
     */
    _SS_ErrorFlag Finalize( MPI_Comm comm /**< communicator for this feature */ );

    /** Return the value calculated by the feature */
    _SS_ErrorFlag Get( double &value /**< output, the return value */ ) const ;

private:

    int status; /**< current status for the feature */
    double final_value /**< the actual value */;

    /**< Initalization function for the feature. */
    virtual _SS_ErrorFlag InitializeFeatureCollection( MPI_Comm comm,
                                                       const int &num_rows /**< number of rows */,
                                                       const int &num_cols /**< number of cols */) = 0;

    /**< The next element function is called in a loop over the non-zeros. */
    virtual _SS_ErrorFlag NextElement( MPI_Comm comm /**< communicator (might not be safe to use FIXME) */ ,
                                       const int &row /**< row index */,
                                       const int &col /**< column index */,
                                       const double &data /**< matrix value */ ) = 0;

    /**< Finalize the feature collection and return the final value. */
    virtual _SS_ErrorFlag FinalizeFeatureCollection(MPI_Comm comm /**< communicator */,
                                                    double &final_value /**< final value */ ) = 0;

};

/** This class represents a collection of features. */
class _SS_Features
{
public:
    MPI_Comm _comm; /**< communicator for this list of features. This is the matrix communicator */
    std::vector< std::shared_ptr<_SS_Feature> > features_list;

    /** Constructor */
    _SS_Features( MPI_Comm comm /**< communicator*/ );

    /** Add a feature to the list  */
    _SS_ErrorFlag AddFeature( std::shared_ptr<_SS_Feature> feature /**< shared_ptr to the feature. */ );

    /**< Call initialize on all features in the list */
    _SS_ErrorFlag Initialize( const int &num_rows, /**< number of rows in the matrix */
                              const int &num_cols  /**< number of columns in the matrix */);

    /**< Call Finalize on all the features in the list */
    _SS_ErrorFlag Finalize();

    /**< Call Next on all the features in the list */
    _SS_ErrorFlag Next(const int &row /**< row index */ , const int &col /**< col index */, const double &value /**< matrix value */);

    /**< Get the final result for all features in the list */
    _SS_ErrorFlag Get( std::map < std::string, double > &fmap /**< map containing the final results */ ) const ;

};


#endif


