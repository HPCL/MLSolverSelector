
#include "CPUTime.h"


using namespace Measurements;

CPUTime::CPUTime( double bad_value ) : _SS_Measurement(bad_value)
{
    parameter_list.insert(std::make_pair("converged",-9999));
}

_SS_ErrorFlag CPUTime::StartMeasurement( MPI_Comm comm, std::map<std::string,double> &mstruct ) 
{
    mstruct["converged"] = -9999;
    time_start = std::chrono::high_resolution_clock::now();
    return _SS_error_flag;
}

_SS_ErrorFlag CPUTime::StopMeasurement(MPI_Comm comm, std::map<std::string,double> &mstruct , double &final_value )
{
    double converged = mstruct["converged"];
    std::cout << converged << " is reason converged \n ";
    if (converged > 0 )
    {
        auto now = std::chrono::high_resolution_clock::now();
        auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(now-time_start).count();
        final_value = (double) durr;
    }
    else
    {
        final_value = std::numeric_limits<double>::max();
    }
    return _SS_error_flag;
}
