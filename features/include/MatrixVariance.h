#ifndef __SS_VARIANCE_HEADER
#define __SS_VARIANCE_HEADER

#include "Features.h"

/* Example feature for calculating the varience */
namespace Features
{

class Variance : public _SS_Feature
{
public:
    /* variance is 1/m ( sum j = 1 : m ) ( a(i,j) - u ) ^ 2 where u is 1/m * the row sum */
    /* to do in one loop, need to expand 1/m (sum) a(i,j)^2 - 2*u + u^2
     * so keep track of a^2 and u */
    double m;
    std::vector<double> squared_sum, sum;
    unsigned int _type;

    Variance( const unsigned int  &type, const std::string &name );
private:
    _SS_ErrorFlag InitializeFeatureCollection( MPI_Comm comm, const int &num_rows, const int &num_cols ) override;

    _SS_ErrorFlag FinalizeFeatureCollection( MPI_Comm comm, double &final_value ) override;

    _SS_ErrorFlag NextElement( MPI_Comm comm, const int &rindex, const int &cindex, const double &value ) override;

};
}
#endif

