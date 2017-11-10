#include "Measurements.h"


/** \file Measurements.C
 * \brief Base classes for measurements.
 **/

_SS_Measurement::_SS_Measurement( )
{
    bad = -1;
    status = 0;
}

_SS_Measurement::_SS_Measurement( double bad_value /**< the bad value for this measurment */ )
{
    bad = bad_value;
    status = 0; // 0 initialized, 1 started, 2 finalized
}

_SS_ErrorFlag _SS_Measurement::Start( MPI_Comm comm /**< communicator*/,
                                      std::map<std::string,double> &mstruct /**< data structure */)
{
    if ( status == 0 || status == 2 )
    {
        final_value = -1;
        StartMeasurement(comm , mstruct );
        status = 1;
    }
    else
        std::cout << " Measurement in wrong state " << status << " \n";

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurement::Stop( MPI_Comm comm /**< communicator */ ,
                                     std::map<std::string, double> &mstruct /**< data struct containing solver data (hopefully) */ )
{
    if ( status == 1 )
    {
        StopMeasurement(comm, mstruct,  final_value );
        status = 2;
    }
    else
        std::cout << " Stop measurement called in status " << status << "\n";

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurement::Get( double &value /**< output, the final value */ ) const
{
    if ( status == 2 )
    {
        value = final_value;
    }
    else
        std::cout << " GetFinalValue called in status " << status << "\n";
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurement::GetBad( double &_bad_value /**< output, the bad value */ ) const
{
    _bad_value = bad;
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurement::GetParameterList( std::map<std::string,double> &list /**< output, list of parameters needed by this solver */ )
{
    list = parameter_list;
    return _SS_error_flag;
}


_SS_Measurements::_SS_Measurements( MPI_Comm comm  /**< communicaator */)
{
    _comm = comm;
}

_SS_ErrorFlag _SS_Measurements::AddMeasurement( const std::string name, std::shared_ptr<_SS_Measurement> measurement )
{
    measurements_list.insert(std::make_pair( name, measurement ) );
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::RemoveMeasurement( const std::string name )
{
    measurements_list.erase( name );
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::Start( std::map< std::string, double > &mstruct /**< the solver data */ )
{
    for (auto &it : measurements_list )
    {
        it.second->Start(_comm, mstruct);
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::Stop( std::map< std::string, double > &mstruct /**< the data structure */  )
{
    for (auto &it : measurements_list )
    {
        it.second->Stop(_comm, mstruct );
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::Get( std::map < std::string, double > &mmap /**< ouput, map to values */) const
{
    double val;
    for (auto &it : measurements_list )
    {
        it.second->Get( val );
        mmap.insert( std::make_pair( it.first, val ) );
    }

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::GetBad( std::map < std::string, double > &mmap /**< output, the bad values */) const
{
    double val;
    for (auto &it : measurements_list )
    {
        it.second->GetBad( val );
        mmap.insert( std::make_pair( it.first, val ) );
    }

    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Measurements::GetParameterList( std::map< std::string, double > &mlist /**< output, blah */ ) const
{
    std::map<std::string, double > mml;
    for ( auto &it : measurements_list )
    {
        mml.clear();
        it.second->GetParameterList( mml );
        for ( auto itt : mml )
            mlist.insert( itt );
    }
    return _SS_error_flag;
}

