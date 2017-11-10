#include "MatrixVariance.h"

/* Example feature for calculating the varience */
using namespace Features;

  Variance::Variance( const unsigned int &type, const std::string &name ) : _SS_Feature(name)
  { 
    _type = type;
  }
  
_SS_ErrorFlag Variance::InitializeFeatureCollection( MPI_Comm comm, const int &num_rows, const int &num_cols ) 
  {

      if ( _type != 1)
      {
        m = num_cols;
        squared_sum.resize(num_rows,0); /* want max over all rows */
        sum.resize(num_rows,0);
      }
      else
      {
        m = num_rows;
        squared_sum.resize(num_cols,0);
        sum.resize(num_cols,0);
      }      
      return _SS_error_flag;
  }

  _SS_ErrorFlag Variance::FinalizeFeatureCollection( MPI_Comm comm, double &final_value ) 
  {
     double glo_max = -1;
     for ( unsigned int i = 0; i < sum.size(); i++ )
     {  
        squared_sum[i] = (1.0/m)*squared_sum[i] - sum[i]*sum[i]/(m*m);
     }
     double loc_max = *(std::max_element(squared_sum.begin(),squared_sum.end()));
     MPI_Allreduce(&loc_max, &glo_max, 1, MPI_DOUBLE, MPI_MAX,comm);
     final_value = glo_max;
     return _SS_error_flag;
  }

  _SS_ErrorFlag Variance::NextElement( MPI_Comm comm, const int &rindex, const int &cindex, const double &value ) 
  {
     int index= ( _type == 1 ) ? rindex : cindex ;
     squared_sum[index] += value*value;
     sum[index] += value;
     return _SS_error_flag;
  }

