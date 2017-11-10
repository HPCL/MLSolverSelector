#ifndef _SS_CPUTIME_HEADER
#define _SS_CPUTIME_HEADER

#include "Measurements.h"

/** An example measurement calculating the CPUTime taken during the solve. TODO allow the user
 * to set the order in which measurements are calculated ! */
namespace Measurements
{

class CPUTime : public _SS_Measurement
{
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> time_start;

    /**Initialize the measurement. In this case, we need to know if the measurement
     * convered, so we add a request for "converged" to the data struct. */
    CPUTime( double bad_value );

private:

    /** Start the CPU Time measurement. Here we set the value of converged in the data struct
     * to -9999. That way, we can make sure it was set by the user. If its not set, we can't be
     * sure what actually happened. Finally, the start time is set. */
    _SS_ErrorFlag StartMeasurement( MPI_Comm comm, std::map<std::string,double> &mstruct ) override;

    /** Stop the measurment. First, we check that the user set the converged parameter, and that it did in fact
     * converge. If it did not, the solve time is set to DOUBLE_MAX, indicating this is an awful solver. If it did,
     * we stop the clock and set the final value. */
    _SS_ErrorFlag StopMeasurement(MPI_Comm comm, std::map<std::string,double> &mstruct , double &final_value ) override;
};
}

#endif
