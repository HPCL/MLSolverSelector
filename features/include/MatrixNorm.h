#ifndef _SS_NORM_HEADER
#define _SS_NORM_HEADER
#include "Features.h"
namespace Features
{

/** An example feature, in this case, the norm */
class Norm : public _SS_Feature
{
public:
    unsigned int _type;
    std::vector< double > data;

    Norm( const unsigned int &type, const std::string &name);

private:
    _SS_ErrorFlag InitializeFeatureCollection(MPI_Comm comm, const int &num_rows, const int &num_cols ) override;
    _SS_ErrorFlag FinalizeFeatureCollection( MPI_Comm comm, double &final_value ) override;
    _SS_ErrorFlag NextElement( MPI_Comm comm, const int &rindex, const int &cindex, const double &value ) override;
};

}
#endif

