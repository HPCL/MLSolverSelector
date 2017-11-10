#ifndef _SS_MEASURE_H
#define _SS_MEASURE_H

#include "ParameterBase.h"


/** \file SolverSelecter.I
 * \brief Base classes for measurements.
 **/


/** Base class for a Measurement */
class _SS_Measurement
{
private:

    int status; /**< status of the measurement */
    double final_value; /**< final value of the measurement */
    double bad; /**< this is the "bad" multiplier for classification of the solvers based on this measurement */;

    /** Virtual function called prior to a solve */
    virtual _SS_ErrorFlag StartMeasurement( MPI_Comm comm /**<communicator */,
                                            std::map<std::string,double> &mstruct /**< requested solver info */ ) = 0;

    /** virtual function called right after the solve. TODO, need to use an ordered list so the measurements
     * are run in order. That way we can start timers last and stop them first. */
    virtual _SS_ErrorFlag StopMeasurement(MPI_Comm comm, /**< communicator */
                                          std::map< std::string, double> &mstruct , /**<requested solver info */
                                          double &value  /**< the final value */) = 0;

public:
    std::map<std::string,double> parameter_list /**< list of the parameters that this measurement needs to know to calculate the result.
  These must be set during the solve by the user. One example is that we need to know if the solver converged in order to assign
  a time. If the solver did not converge, the solve time is useless and should be set as such.  */ ;

    _SS_Measurement( ) ;
    _SS_Measurement( double bad_value /**< the bad value for this measurment */ );

    /** Check the status and start the measurement. The measurement should add any required parameters
     * to the mstruct.  */
    _SS_ErrorFlag Start( MPI_Comm comm /**< communicator*/,
                         std::map<std::string,double> &mstruct /**< data structure */);

    /** Check the status and stop the measurement. Mstruct is a data structure set by
     * the user during the solve. This contains all the data needed by the measurements. For
     * example, this might contain information on if the solver converged.  */
    _SS_ErrorFlag Stop( MPI_Comm comm /**< communicator */ ,
                        std::map<std::string, double> &mstruct /**< data struct containing solver data (hopefully) */ );

    /** Get the final value for the measurement */
    _SS_ErrorFlag Get( double &value /**< output, the final value */ ) const ;

    /** Get the bad multiplier for this measurement */
    _SS_ErrorFlag GetBad( double &_bad_value /**< output, the bad value */ ) const ;

    /** get the parameter list for this measurement */
    _SS_ErrorFlag GetParameterList( std::map<std::string,double> &list /**< output, list of parameters needed by this solver */ );


};

/** This class represents a collection of measurements. */
class _SS_Measurements
{
public:

    MPI_Comm _comm /**< MPI communicator */ ;
    std::map< std::string, std::shared_ptr<_SS_Measurement> > measurements_list; /**< list of measurements */

    /** Constructor */
    _SS_Measurements( MPI_Comm comm  /**< communicaator */) ;

    /** Add a new measurement to the list. */
    _SS_ErrorFlag
    AddMeasurement( const std::string name, std::shared_ptr<_SS_Measurement> measurement );

    /**< Remove a measurement from the list */
    _SS_ErrorFlag
    RemoveMeasurement( const std::string name );

    /** Call the StartMeasuring function on all measurements in the list */
    _SS_ErrorFlag
    Start( std::map< std::string, double > &mstruct /**< the solver data */ );

    /** Call the stop measuring function on all the measurements in the list */
    _SS_ErrorFlag
    Stop( std::map< std::string, double > &mstruct /**< the data structure */  );

    /**< Get the values calculated by the measurements in the list */
    _SS_ErrorFlag
    Get( std::map < std::string, double > &mmap /**< ouput, map to values */) const ;

    /**< Get a map of the bad values for every measurement in the list */
    _SS_ErrorFlag
    GetBad( std::map < std::string, double > &mmap /**< output, the bad values */) const ;

    /** Get a map of the parameter lists required by each measurement */
    _SS_ErrorFlag
    GetParameterList( std::map< std::string, double > &mlist /**< output, blah */ ) const ;
};



#endif

