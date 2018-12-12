
#ifndef PETSCFE_HEADER
#define PETSCFE_HEADER

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include <memory>
#include <mpi.h> // MPI Support
#include <algorithm>
#include "petscdmda.h"
#include "petscksp.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

// This will extract features from any PETSC matrix using mat-vecs and the sampling system. All features 
// are seperated out using MACROS, so to change the feature set, we only need to change the macro definition
// and recompile. This gets rid of a tonne of run-time if statements. 
//
// Symmetry based Features: Symmetry is somewhat finicky to implement in parallel. To avoid that, I implemented the symmetrical
// features by gathering the sample data using a MPI sum inside the above all reduce command. Because we only
// use sample data for these calcualtions, this equates to sending NUMSAMPLES**2 extra doubles in the communication buffer.
// I don't see this ever being a problem unless the number of samples grows dramatically, in which case we will have other 
// problems ( i.e, the mat-vecs ) . 


// Define a tolerance for a value being non-zero
#define NONZEROTOLERANCE 1e-125

// Define a processor to do the final calcuations on 
#define ROOTPROC 0

// Defines for features supported. The value associated with each feature is the number of entries
// used in the MPI buffer to do the communication. I.e, diagonal sign needs three places in the buffer
// to make sure we can calculate the value of the root processor. 

// MPI_SUM
#define _NNZ 1 
#define _AVGNONZEROSPERROW 1
#define _ABSOLUTENONZEROSUM 1
#define _FROBENIUSNORM 1 
#define _AVERAGEDIAGONALDISTANCE 2                 
#define _ROWVARIANCE 1
#define _TRACE 1
#define _ABSOLUTETRACE 1
#define _DIAGONALMEAN 1          
#define _DIAGONALAVERAGE 2    
#define _DIAGONALNONZEROS 1 

//MPI_MAX
#define _INFINITYNORM 1
#define _MAXNONZEROSPERROW 1
#define _DIAGONALSIGN 3           
#define _LOWERBANDWIDTH 1 
#define _UPPERBANDWIDTH 1 
#define _ROWDIAGONALDOMINANCE 1

//MPI_MIN
#define _MINNONZEROSPERROW 1

//MPI_COLUMN
#define _ONENORM 1
#define _COLDIAGONALDOMINANCE 1
#define _COLUMNVARIANCE 2 

//MPI_SYMMETRY
#define _SYMMETRICITY 1
#define _NONZEROPATTERNSYMMETRY 1
#define _SYMMETRICINFINITYNORM 1
#define _SYMMETRICFROBENIUSNORM 1
#define _ANTISYMMETRICINFINITYNORM 1
#define _ANTISYMMETRICFROBENIUSNORM 1



#include "FeatureSets.h"


// Calculate the total number of data points needed in comm vector. 
#define FEATURECOUNT(npoints) MAXCOUNT + MINCOUNT + COLCOUNT*npoints + SUMCOUNT + SYMMETRY*npoints*npoints 

int MapDMNumberingToNaturalNumbering( std::vector< std::pair< int , int > > &points ) 
{
    // TODO Need to make sure we are still sampleing middle rows and edge rows. Other than 
    // that I think everything else just follows through. 
    return 0;
    
}

// TODO This function uses set to avoid a custom sort routine, then copies that set over into
// a vector. Either need to impliment custom sort, or change code to use set instead of vector. 
int GetSamplePoints( int n, int edge, int interior , std::vector<std::pair<int,int>> &points_vec )
{
    srand(time(NULL));
    
    std::set<std::pair<int,int>> points; 
    if ( edge < 0 && interior < 0 ) {
       // Do all interior points 
       for ( int i = 0; i < n ; i++ ) {
         points.insert( std::make_pair( i, 1 ) );
      } 
  
    }
    else if ( edge + interior >= n ) 
    {
      // Do all rows with some edges
      for ( int i = 0; i < n ; i ++ ) 
       points.insert( std::make_pair(i, ( i > n - 3 || i < 3 ) ? -1 : 1 ) ) ;
    }
    else 
    {     
        int edge_min, edge_max;

        for ( int i = 0; i < edge / 2 ; i++ ) {
           points.insert( std::make_pair( i, -1 ) );
           points.insert( std::make_pair( n-i-1, -1 ) );
        }

        if ( edge % 2 == 1 ) { 
          points.insert( std::make_pair( n - (edge)/2 -1 , -1 ) );    
        }
        
        edge_max = ( edge % 2 == 1 ) ? n - (edge)/2 - 1 : n - edge/2 - 2 ;
        edge_min = edge/2 -1 ;
        int failsafe = 0;
        
        //TODO -- This is junk. We cycle through random numbers until we find 
        //one that isnt in the vector aleady. 
        while ( points.size() < interior + edge && failsafe++ <= 100*n  ) {
            int ind = rand() % n ;
            if ( ind < edge_max && ind > edge_min ) { 
              points.insert(std::make_pair(ind, 1) );
            }
        }

    }

    points_vec.clear();
    points_vec.assign( points.begin(), points.end()); 

    return 0; 
}

int GetJacobianColumn( Mat J, std::pair< int, int > &point, Vec *s, bool mat_vec=true )
{

    if (!mat_vec)
    {
      Vec basis, extra;
      MatCreateVecs(J,&basis, &extra);
      VecDuplicate(extra, s);
      MatGetColumnVector(J,*s,point.first);
      VecDestroy(&basis);
      VecDestroy(&extra);
    }
    else 
    {    

      // basis is vector
      Vec basis, extra; // The basis vector we will multiply by
      MatCreateVecs( J, &basis, &extra ) ; // make basis compatible with J. 
      VecDuplicate(extra, s);

      PetscInt high, low;  
      VecGetOwnershipRange(basis, &low, &high );
      
      VecSet(basis,0);     
      if ( point.first >= low && point.first < high );
          VecSetValue( basis, point.first, (PetscScalar) 1.0 , INSERT_VALUES ); 
      
      VecAssemblyBegin(basis);
      VecAssemblyEnd(basis);
      VecAssemblyBegin(*s);
      VecAssemblyEnd(*s);
      MatMult( J, basis, *s );  

 
    }

    
    return 0;
}


int GetJacobianSamplePointsPercentage( int n , double percentage , int edge_rows, std::vector<std::pair<int,int>> &points ) {

  // Gets a fixed percentage of the interior rows  ( total - edge_rows ) with a minimum of 1 


  if ( percentage < 1.0 && percentage > 0.0 ) {
      
      int int_rows = ( ( n - edge_rows ) * percentage ) ;

      GetSamplePoints( n, edge_rows, ( int_rows < 1 ) ? 1 : int_rows , points );

  } else {
    std::cout << " \t Percentage must be a fraction between 0 and 1 \n " ;
  }
  return 0;

}

int GetJacobianSamplePoints( Mat J , int edge, int interior, std::vector< std::pair< int, int > >&points )
{
    PetscInt n,m;
    MatGetSize(J, &m, &n ) ;
    if ( edge < 0 && interior < 0 )
      GetSamplePoints( n , edge, interior, points );
    else if ( edge < 0 || interior < 0 )
      GetSamplePoints( n, n, n, points ) ; 
    else  
      GetSamplePoints( n, edge, interior,points ) ; 
}

void MPI_FeatureReduce( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
    double *in = (double*) invec;
    double *inout = (double*) inoutvec; 
   
    // The first SUMCOUNT indicies are just MPISUMS. 
    int shift = 0;
    for ( int i = shift ; i < shift + SUMCOUNT ; i++ )  inout[i] += in[i] ;
        
    // Next are the Maximum features    
    shift = SUMCOUNT;
    for ( int i = shift; i < shift + MAXCOUNT ; i++ )   inout[i] = ( inout[i] > in[i] ) ? inout[i] : in[i] ; 
    
    //Next are the minimums. 
    shift += MAXCOUNT;
    for ( int i = shift; i < shift + MINCOUNT;  i++ )   inout[i] = ( inout[i] < in[i] ) ? inout[i] : in[i] ;
    
    //Everything else is a MPISUM 
    shift += MINCOUNT;
    for ( int i = shift ; i < *len; i++)   inout[i] += in[i] ; 
}

int ExtractJacobianFeatures( Mat J , int edge, int interior, std::map<  std::string, double > &fnames , bool matvecs=true )
{
  std::vector< std::pair< int, int >  >points ; 
  std::map< std::string, double > feature_set; 
  
  PetscInt n, m, low, high;
  MatGetSize(J, &n, &m ) ; 
 
  
  Vec *s ; 
  const PetscScalar *array; 
  
  double nsamples_mid = 0;

  bool sample;
  double val, aval, t,v;
  int type, k;
   
  std::vector< std::pair< int, int > > sample_points;
  GetJacobianSamplePoints( J, edge, interior, points ); 
  int npoints = points.size();


  #ifdef NNZ 
  double nnz = 0;   
  #endif   
 
  #ifdef AVGNONZEROSPERROW 
  double avg_nnz = 0;
  #endif

  #ifdef ABSOLUTENONZEROSUM
  double abs_nz_sum = 0.0; // Sum over all non zeros in the middle block  
  #endif

  #ifdef TRACE   
  double trace = 0.0;     // trace 
  #endif   

  #ifdef ABSOLUTETRACE   
  double abs_trace = 0.0; // absolute value of the trace; 
  #endif

  #ifdef FROBENIUSNORM   
  double fro_norm  = 0.0 ; // all middle points  
  #endif 

  #ifdef DIAGONALAVERAGE  
  double diag_ave = 0.0; // average of abs value of diagonal 
  double diag_ave_denom = 0; // num nz on diag to divide by
  #endif

  #ifdef AVERAGEDIAGONALDISTANCE   
  double ave_diag_dist = 0.0; // average distance of non-zeros from the diagonal -- middle points all  
  double ave_diag_dist_denom = 0 ; 
  #endif

  #ifdef LOWERBANDWIDTH   
  double blower = 0;
  #endif

  #ifdef UPPERBANDWIDTH   
  double bupper = 0;
  #endif

  #if defined(LOWERBANDWIDTH) || defined(UPPERBANDWIDTH) 
  int bwidth = 0;
  #endif 

  #ifdef DIAGONALMEAN 
  double diag_mean = 0.0;
 #endif

  #ifdef DIAGONALSIGN
  double diag_sign = 3; 
  #endif

  #ifdef DIAGONALNONZEROS
  double diag_nnz = 0;  
  #endif

  #ifdef ROWVARIANCE 
  double row_var = 0.0;
  #endif    

  #ifdef MAXNONZEROSPERROW   
  std::vector< int > max_nnz_row(npoints); // nnz_max_rows uses sample points only from the middle only
  #endif

  #ifdef MINNONZEROSPERROW
  //Set all these values to n+1 to start with. Since this is impossible ( i.e, more nz than rows ), we can use that to test weather
  //a row was accessed on this processor. 
  std::vector< int > min_nnz_row(npoints, n+1); // nnz_max_rows uses sample points only from edges only 
  #endif

  #ifdef ONENORM   
  std::vector<double> one_norm(npoints); // This is the abs sum of column for all samples this one needs mpi sum then max 
  #endif

  #ifdef INFINITYNORM   
  std::vector<double> inf_norm(npoints); // All sample points , mid and edges
  #endif

  #ifdef COLDIAGONALDOMINANCE   
  std::vector<double> col_diag_dom(npoints);
  #endif 

  #ifdef ROWDIAGONALDOMINANCE   
  // Row diagonal dominance uses a MPI over the definition for each row, with diagonal dominance only
  // true if every row is negative. So, to make sure this works in parallel, make every row negative,
  //  then set the sum to zero when this row is first accessed.  
  std::vector<double> row_diag_dom(npoints,-1.0);
  #endif   

  #ifdef ROWVARIANCE
  std::vector<double> row_squared_sum(npoints);  // THese are for the variences s^2 _i = ( E(x^2) - n ( E x )^2 )/(n-1) 
  std::vector<double> row_sum(npoints);
  #endif

  #ifdef COLUMNVARIANCE   
  std::vector<double> col_squared_sum(npoints);
  std::vector<double> col_sum(npoints); 
  #endif

  #ifdef DIAGONALSIGN
  int neg(0), pos(0),zero(0); 
  #endif

  #if SYMMETRY    
  std::vector<double> sample_data(npoints * npoints);  // all zeros, set row*npoints + column 
  #endif 
 
   int num_t1 = 0;
   int num_samples = 0;
   for ( int i = 0; i < npoints; i++ ) 
   {
       
       Vec vec;
       GetJacobianColumn( J, points[i], &vec, matvecs );
       
       VecGetArrayRead( vec, &array );
       VecGetOwnershipRange( vec, &low, &high ) ;  
       type = points[i].second;
       if ( type == 1 ) nsamples_mid++;

       k = 0;
       while ( points[k].first < low && k < points.size() ) k++;   // Find the first sample point. 
       for ( int j = low; j < high; j++ ) 
       { 
             val = array[j-low];
             
             aval = abs(val);
                     

              // Is this a Sample point ? 
             if ( k < points.size() && j == points[k].first ) 
             {   
               sample=true; k++;
               #ifdef MINNONZEROSPERROW 
               if (min_nnz_row[k-1] == n+1) {
                  min_nnz_row[k-1] = 0; 
               }
               #endif

               #ifdef ROWDIAGONALDOMINANCE
               row_diag_dom[k-1] = 0.0;
               #endif
               num_samples++; 
             }
             else 
               sample = false ; 
             

             // If this is a middle point 
             if ( type == 1 ) 
             {

                 num_t1 ++; 
                 if ( aval > NONZEROTOLERANCE ) 
                 {
                      
                      #ifdef NNZ     
                      nnz++;
                      #endif
                     
                      #ifdef AVGNONZEROSPERROW  
                      avg_nnz++;
                      #endif                      

                      #ifdef ABSOLUTENONZEROSUM                      
                      abs_nz_sum += aval; 
                      #endif

                      #ifdef FROBENIUSNORM 
                      fro_norm += aval*aval;
                      #endif   

                      #ifdef AVERAGEDIAGONALDISTANCE
                      if ( points[i].first != j )
                      {
                         ave_diag_dist += abs( points[i].first - j ) ; 
                         ave_diag_dist_denom ++; 
                      }
                      #endif

                      #ifdef COLUMNVARIANCE 
                      col_squared_sum[i] += aval*aval;
                      col_sum[i] += val; 
                      #endif
                      
                      if (sample) {
                        
                        #ifdef ROWVARIANCE 
                        row_squared_sum[k-1] += aval*aval;
                        row_sum[k-1] += val; 
                        #endif

                        #ifdef MAXNONZEROSPERROW 
                        max_nnz_row[k-1]++;
                        #endif
                      }
                 }
            }
               
            #ifdef MINNONZEROSPERROW // Count non-zeros in the sample rows. 
            if (aval > NONZEROTOLERANCE && sample ) min_nnz_row[k-1]++;
            #endif
             
             // These are features that use both mid and edge points 
             #ifdef ONENORM
             one_norm[i] += aval;
             #endif

             if (points[i].first == j ) 
             {       
               #ifdef DIAGONALSIGN 
               if ( val < -1.0*NONZEROTOLERANCE ) neg = 1;
               else if ( aval < NONZEROTOLERANCE ) zero = 1;
               else pos = 1;
               #endif

               if ( sample )
               { 
                  #ifdef TRACE
                  trace += val;
                  #endif

                  #ifdef DIAGONALMEAN  
                  diag_mean += val;
                  #endif
                    
                  #ifdef ABSOLUTETRACE
                  abs_trace += aval;
                  #endif
               
                  #ifdef DIAGONALNONZEROS         
                  if ( aval > NONZEROTOLERANCE ) diag_nnz++;
                  #endif
                  
                  // Todo, this is redundent if diag non zeros is defined 
                  #ifdef DIAGONALAVERAGE 
                  if (aval > NONZEROTOLERANCE ) {
                     diag_ave += aval; 
                     diag_ave_denom++;
                  }
                  #endif 

                } 
             } 

             #ifdef COLDIAGONALDOMINANCE             
             col_diag_dom[i] += ( ( points[i].first == j ) ? -1 : 1 )*aval; 
             #endif

             if (sample) {
                   #ifdef INFINITYNORM 
                   inf_norm[k-1] += aval; // points[k-1] = j, so save it in k-1. 
                   #endif

                   #ifdef ROWDIAGONALDOMINANCE                   
                   row_diag_dom[k-1] += ( ( points[i].first == j ) ? -1 : 1 )*aval ;                     
                   #endif

                   #if SYMMETRY                   
                   sample_data[(k-1)*npoints + i ] = val; 
                   #endif             
             }
             
             #if defined(LOWERBANDWIDTH) || defined(UPPERBANDWIDTH)             
             bwidth = points[i].first - j ;
             if ( aval > NONZEROTOLERANCE )
             {
                  #ifdef UPPERBANDWIDTH             
                  if ( bwidth  > bupper ) bupper = bwidth ;
                  #endif
                  #ifdef LOWERBANDWIDTH             
                  if ( -bwidth > blower ) blower = -bwidth ; 
                  #endif
             }
             #endif // LOWERBANDWIDTH || UPPERBANDWIDTH
       }          
       VecRestoreArrayRead( vec , &array ) ;
       VecDestroy(&vec);
   }
   #ifdef ROWVARIANCE   
   row_var = 0.0;
   for ( auto i = 0; i < npoints; i++ ) 
       row_var += ( 1./ ( nsamples_mid*(nsamples_mid-1) ) )* ( row_squared_sum[i] - (1.0/nsamples_mid) * row_sum[i] * row_sum[i] ); 
   #endif 

   // Finish it 
   int feature_count  = FEATURECOUNT(npoints); 
   
   
   double *features = (double*) malloc( sizeof(double) * ( feature_count  ) );
   double *rfeatures = (double*) malloc( sizeof(double) * ( feature_count ) ) ;
   int c = 0;

   #ifdef NNZ   
   features[c++] =  nnz * n / (double) nsamples_mid ;                  
   #endif

   #ifdef AVGNONZEROSPERROW   
   features[c++] = avg_nnz / (double) nsamples_mid ;                       
   #endif

   #ifdef ABSOLUTENONZEROSUM   
   features[c++] = abs_nz_sum * n / (double) nsamples_mid;             
   #endif

   #ifdef TRACE   
   features[c++] = trace * n / (double) npoints  ;
   #endif
   
   #ifdef ABSOLUTETRACE   
   features[c++] = abs_trace * n / (double) npoints  ;
   #endif
 
   #ifdef DIAGONALMEAN  
   features[c++] = diag_mean / (double) npoints ;
   #endif
   
   #ifdef DIAGONALNONZEROS 
   features[c++] = diag_nnz * n / (double) npoints ;  
   #endif

   #ifdef DIAGONALAVERAGE   
   features[c++] = diag_ave;
   features[c++] = diag_ave_denom;   
   #endif 
   
   #ifdef FROBENIUSNORM   
   features[c++] = fro_norm * n * n / ( nsamples_mid * nsamples_mid ) ;
   #endif 

   #ifdef ROWVARIANCE    
   features[c++] = row_var;               
   #endif 

   #ifdef AVERAGEDIAGONALDISTANCE
   features[c++] = ave_diag_dist ;
   features[c++] = ave_diag_dist_denom ;                                
   #endif

   /*  Start the MPI MAX FEATURES. *****************************************/  
   #ifdef INFINITYNORM   
   features[c++] = *std::max_element(inf_norm.begin(),inf_norm.end());  
   #endif
   
   #ifdef MAXNONZEROSPERROW   
   features[c++] = *max_element(max_nnz_row.begin(), max_nnz_row.end());
   #endif

   #ifdef DIAGONALSIGN   
   features[c++] = neg;                                                 
   features[c++] = pos;
   features[c++] = zero;
   #endif

   #ifdef LOWERBANDWIDTH   
   features[c++] = blower;                                              
   #endif

   #ifdef UPPERBANDWIDTH   
   features[c++] = bupper;                                              
   #endif

   #ifdef ROWDIAGONALDOMINANCE  
   features[c++] = *std::max_element(row_diag_dom.begin(), row_diag_dom.end()); 
   #endif

   /* Start the MPI MIN Features ********************************************/ 
   #ifdef MINNONZEROSPERROW  
   features[c++] = *min_element(min_nnz_row.begin(), min_nnz_row.end());   
   #endif 

   /* Start the column sum features ******************************************/
   #ifdef ONENORM 
   for ( auto &it : one_norm ) features[c++] = it;
   #endif

   #ifdef COLDIAGONALDOMINANCE   
   for ( auto &it : col_diag_dom ) features[c++] = it;
   #endif
   
   #ifdef COLUMNVARIANCE   
   for ( auto &it : col_squared_sum ) features[c++] = it;
   for ( auto &it : col_sum ) features[c++] = it;
   #endif

   /* Add the sample data ***************************************************/ 
   #if SYMMETRY   
   for ( auto &it : sample_data ) features[c++] = it;
   #endif 



   // MPI Reduction -- We actually want this to be a simple Reduce, right? It would be catastrophic
   // if the Solver selected on two processors was different, so we will have to do a check anyway. May
   // as well just do these calculations and classification on the root then scatter the solver choice. (Unless
   // ML models are always deterministic ? ) 
   
   
   MPI_Op custom_mpi_operation;
   MPI_Op_create( MPI_FeatureReduce, true, &custom_mpi_operation );
   MPI_Reduce( features, rfeatures, feature_count, MPI_DOUBLE, custom_mpi_operation, ROOTPROC, PETSC_COMM_WORLD );

   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
   if ( rank == ROOTPROC )
   { 
     //Write the feature set. 
     fnames.clear();

     fnames.insert( std::make_pair( "Dimension",             n  ) );
     fnames.insert( std::make_pair( "npoints",               npoints  ) );
     
     c = 0;    
     #ifdef NNZ 
     fnames.insert( std::make_pair( "nnz",                  rfeatures[c++] ) );
     #endif   

     #ifdef AVGNONZEROSPERROW 
     fnames.insert( std::make_pair( "AvgNonzerosPerRow",    rfeatures[c++]  ) );
     #endif

     #ifdef ABSOLUTENONZEROSUM
     fnames.insert( std::make_pair( "AbsoluteNonZeroSum",   rfeatures[c++]  ) );
     #endif   

     #ifdef TRACE
     fnames.insert( std::make_pair( "Trace",                rfeatures[c++]  ) );
     #endif

     #ifdef ABSOLUTETRACE 
     fnames.insert( std::make_pair( "AbsoluteTrace",        rfeatures[c++]  ) );
     #endif   
     
     #ifdef DIAGONALMEAN    
     fnames.insert( std::make_pair( "DiagonalMean",         rfeatures[c++]  ) );
     #endif

     #ifdef DIAGONALNONZEROS 
     fnames.insert( std::make_pair( "DiagonalNonZeros",     rfeatures[c++]  ) );
     #endif

     #ifdef DIAGONALAVERAGE  
     fnames.insert( std::make_pair( "DiagonalAverage",      rfeatures[c++] / rfeatures[c++] ) );
     #endif
   
     #ifdef FROBENIUSNORM   
     fnames.insert( std::make_pair( "FrobeniusNorn",        sqrt(rfeatures[c++])) );
     #endif  
     
     #ifdef ROWVARIANCE
     fnames.insert( std::make_pair( "RowVariance",          rfeatures[c++] ) );
     #endif
     
     #ifdef AVERAGEDIAGONALDISTANCE   
     fnames.insert( std::make_pair( "AvgDiagonalDistance" , rfeatures[c++]/ rfeatures[c++] ) );
     #endif
     
     #ifdef INFINITYNORM
     fnames.insert( std::make_pair( "InfinityNorm",         rfeatures[c++] ) );
     #endif
     
     #ifdef MAXNONZEROSPERROW
     fnames.insert( std::make_pair( "MaxNonZeroPerRow",     rfeatures[c++] ) );
     #endif   

     #ifdef DIAGONALSIGN
     // Perform all the buffer post-processing .. 
     neg = rfeatures[c++]; pos = rfeatures[c++]; zero = rfeatures[c++];   
     if ( neg && pos ) diag_sign = -3;
     else if ( neg && zero ) diag_sign = -1; 
     else if ( pos && zero ) diag_sign = -1;
     else if ( neg ) diag_sign = -2; 
     else if ( pos ) diag_sign = 2;
     else if ( zero) diag_sign = 0;  
     fnames.insert( std::make_pair( "DiagonalSign",         diag_sign ) );
     #endif 

     #ifdef LOWERBANDWIDTH   
     fnames.insert( std::make_pair( "lowerBandwidth",       rfeatures[c++] ) );
     #endif

     #ifdef UPPERBANDWIDTH
     fnames.insert( std::make_pair( "upperBandwidth",       rfeatures[c++] ) );
     #endif
     
     #ifdef ROWDIAGONALDOMINANCE
     fnames.insert( std::make_pair( "RowDiagonalDominance", ( rfeatures[c++] < 0 ) ? 1 : 0 ) );
     #endif

     #ifdef MINNONZEROSPERROW
     fnames.insert( std::make_pair( "MinNonZeroPerRow",     rfeatures[c++] ) );
     #endif

     #ifdef ONENORM
     double fone_norm = -1;
     for ( int i = 0; i < npoints; i++ ) {
       t = rfeatures[c++];
       if ( t > fone_norm ) fone_norm = t;
     }
     fnames.insert( std::make_pair( "OneNorm", fone_norm ) );
     #endif

     #ifdef COLDIAGONALDOMINANCE
     double ct; 
     double fccd = 1; 
     for ( int i = 0; i < npoints; i++ ) {
        t = rfeatures[c++];
        if ( t > 0 ) fccd = 0;
     }
     fnames.insert( std::make_pair( "ColumnDiagonalDominance", fccd ) );
     #endif   

      #ifdef COLUMNVARIANCE
      double cvar = 0; 
      for ( int i = 0; i < npoints; i++ ) {
        t = rfeatures[c];  // squared sum over all rows ( n of them ) 
        v = rfeatures[c + npoints];   // sum over all rows ( n of them ) 
        
        // mean = 1/#samples Sum ( var ) 1/npoints  = 
        double ih = ( ( 1./ (n-1) ) )*( t - (1.0/n) * v * v ); 
        cvar += ( 1.0 / (double) nsamples_mid ) * ih; 
        c++;
      }
      c += npoints;
      fnames.insert( std::make_pair( "ColumnVariance", cvar  ) );
      #endif

     #if SYMMETRY
     bool sym = true;
     bool nzsym = true;
     double sfronorm(0),asfronorm(0);
     std::vector<double> sinfnorm(npoints),asinfnorm(npoints);
     for ( int k = 0; k < npoints*npoints; k++ )
     {
          int row = k / npoints;
          int col = k % npoints;
          double value  = rfeatures[ c + row*npoints + col ] ;
          double svalue = rfeatures[ c + col*npoints + row ] ;
          
          #ifdef SYMMETRICITY        
          sym = ( sym && abs( value - svalue ) < 1e-13 );
          #endif

          #ifdef NONZEROPATTERNSYMMETRY        
          nzsym = ( nzsym && (  ( abs(value) > 1e-13 && abs(svalue) > 1e-13 ) || ( abs(value) < 1e-13 && abs(svalue) < 1e-13 ) ) );
          #endif

          #ifdef SYMMETRICINFINITYNORM       
          sinfnorm[row] += abs( 0.5 * ( svalue + value ) );
          #endif
    
          #ifdef ANTIYSYMMETRICINFINITYNORM
          asinfnorm[row] += abs( 0.5 * ( svalue - value ) );
          #endif  

          #if defined(SYMMETRICFROBENIUSNORM) || defined(ANTISYMMETRICFROBENIUSNORM)
          if ( points[col].second == 1 ) {
            #ifdef SYMMETRICFROBENIUSNORM            
            sfronorm += ( 0.5 * ( svalue + value ) ) * ( 0.5 * ( svalue + value ) ) ; 
            #endif

            #ifdef ANTISYMMETRICFROBENIUSNORM  
            asfronorm += ( 0.5 * ( svalue - value ) ) * ( 0.5 * ( svalue - value ) ) ; 
            #endif       
          }
          #endif
     
     }
     #ifdef SYMMETRICITY
     fnames.insert( std::make_pair( "Symmetricity", sym  ) );
     #endif
    
     #ifdef NONZEROPATTERNSYMMETRY
     fnames.insert( std::make_pair( "NonZeroPatternSymmetry", nzsym ) );
     #endif 
    
     #ifdef SYMMETRICINFINITYNORM   
     fnames.insert( std::make_pair( "SymmeticInfinityNorm", *std::max_element(sinfnorm.begin(), sinfnorm.end()) ));
     #endif
    
     #ifdef ANTISYMMETRICINFINITYNORM
     asfronorm = n * sqrt(asfronorm) / (double) nsamples_mid ;
     fnames.insert( std::make_pair( "AntiSymmetricFrobeniusNorm" , asfronorm ) );
     #endif
    
     #ifdef SYMMETRICFROBENIUSNORM   
     sfronorm = n * sqrt(sfronorm) / (double) nsamples_mid ;
     fnames.insert( std::make_pair( "SymmetricFrobeniusNorm", sfronorm   ) );
     #endif
    
     #ifdef ANTISYMMETRICFROBENIUSNORM   
     asfronorm = n * sqrt(asfronorm) / (double) nsamples_mid ;
     fnames.insert( std::make_pair( "AntiSymmerticInfinityNorm", *std::max_element(asinfnorm.begin(),asinfnorm.end())) );
     #endif

     #endif  // SYMMETRY
     
   


   }  

   free(features);
   free(rfeatures);  
   return 0;
}
   
         
#endif 

