#include "Features.h"

/**
 * \file Features.C
 * \brief Impliment file for a "feature".
 **/


_SS_Feature::_SS_Feature(const std::string &_name)
{
    name = _name; 
    status = 0;
    final_value = -1;
}

_SS_ErrorFlag _SS_Feature::Initialize( MPI_Comm comm /**< MPI communicator to be used with this feature */,
                                       const int &num_rows /**< number of rows in the matrix */,
                                       const int &num_cols /**< number of columns in the matrix */)
{
    InitializeFeatureCollection( comm, num_rows, num_cols );
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Feature::Next( MPI_Comm comm /**< communicator for this feature*/,
                                 const int &row /**< row index */,
                                 const int &col /**< col index */,
                                 const double &data /**< data value */ )
{
    NextElement(comm, row, col, data);
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Feature::Finalize( MPI_Comm comm /**< communicator for this feature */ )
{
    FinalizeFeatureCollection( comm, final_value );
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Feature::Get( double &value /**< output, the return value */ ) const
{
    value = final_value;
    return _SS_error_flag;
}


_SS_Features::_SS_Features( MPI_Comm comm /**< communicator*/ )
{
    _comm = comm;
}

_SS_ErrorFlag _SS_Features::AddFeature( std::shared_ptr<_SS_Feature> feature /**< shared_ptr to the feature. */ )
{
    features_list.push_back( feature );
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Features::Initialize( const int &num_rows, /**< number of rows in the matrix */
                                        const int &num_cols  /**< number of columns in the matrix */)
{
    
    for (auto &it : features_list )
    {
        it->Initialize(_comm, num_rows, num_cols);
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Features::Finalize()
{
    for (auto &it : features_list )
    {
        it->Finalize(_comm);
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Features::Next(const int &row /**< row index */ ,
                                 const int &col /**< col index */,
                                 const double &value /**< matrix value */)
{
    for ( auto it : features_list )
    {
        it->Next( _comm, row, col, value );
    }
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_Features::Get( std::map < std::string, double > &fmap /**< map containing the final results */ ) const
{
    double val;
    for ( auto &it :features_list )
    {
        it->Get( val );
        fmap.insert( std::make_pair( it->name, val ) );
    }
    return _SS_error_flag;
}



