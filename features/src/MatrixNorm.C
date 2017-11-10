#include "MatrixNorm.h"

using namespace Features;

Norm::Norm( const unsigned int &type, const std::string &name) : _SS_Feature(name)
{
    _type = type;
}

_SS_ErrorFlag Norm::InitializeFeatureCollection(MPI_Comm comm, const int &num_rows, const int &num_cols ) 
{
    if ( _type == -1)
        data.resize(num_rows,0);
    else if ( _type == 1)
        data.resize(num_cols,0);
    else
        data.resize(1,0);
    return _SS_error_flag;
}

_SS_ErrorFlag Norm::FinalizeFeatureCollection( MPI_Comm comm, double &final_value ) 
{
    if ( _type == -1 || _type == 1 )
    {
        double glo_max = -1;
        double loc_max = *(std::max_element(data.begin(),data.end()));
        //MPI_Allreduce(&loc_max, &glo_max, 1, MPI_DOUBLE, MPI_MAX,comm);
        final_value = loc_max;
    }
    else
    {
        double glo_sum = -1;
        double loc_sum = data[0];
        //MPI_Allreduce(&loc_sum, &glo_sum, 1, MPI_DOUBLE, MPI_SUM, comm );
        final_value = loc_sum ; /* FIXME sqrt */
    }

    return _SS_error_flag;
}

_SS_ErrorFlag Norm::NextElement( MPI_Comm comm, const int &rindex, const int &cindex, const double &value ) 
{
    if ( _type == -1)
        data[rindex] += abs(value);
    else if ( _type == 1)
        data[cindex] += abs(value);
    else
        data[0] += value*value; /*Frob*/

    return _SS_error_flag;
}
