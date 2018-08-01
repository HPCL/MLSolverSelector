
#include "PetscFeatures.h"

#if WITH_PETSCUI 

namespace SolverSelecter
{

PetscTestingSpace::PetscTestingSpace() {}



void
DefaultPetscTestingSpace::extract_features(KSP &ksp,
        std::map<std::string, double> &fmap )
{
    Mat AA,PP;
    KSPGetOperators(ksp, &AA, &PP );
    ExtractJacobianFeatures(AA,10,10,fmap);

    int comm_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &comm_size) ;
    fmap["comm_size"] = (double) comm_size;
}

void
DefaultPetscTestingSpace::start_measurements(KSP &ksp,
        Vec &x,
        Vec &b )

{
    time_start = std::chrono::high_resolution_clock::now();
}

void
DefaultPetscTestingSpace::stop_measurements(KSP &ksp,
        Vec &x,
        Vec &b,
        std::map<std::string, double> &mmap )
{
    auto now = std::chrono::high_resolution_clock::now();
    auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(now-time_start).count();
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    mmap["CPUTime"] =  ( reason > 0 ) ? (double) durr : -1.0 ;
}

void
DefaultPetscTestingSpace::classify_measurements( std::map<std::string, double> &cmap )
{
    cmap["CPUTime"] = 0.3;
}

int
DefaultPetscTestingSpace::MapDMNumberingToNaturalNumbering( std::vector< std::pair< int , int > > &points )
{
    return 0;

}

int
DefaultPetscTestingSpace::GetSamplePoints(int n,
        int edge,
        int interior,
        std::vector<std::pair<int,int>> &points )
{
    srand(time(NULL));
    points.clear();
    if ( edge < 0 && interior < 0 )
    {
        // Do all interior points
        for ( int i = 0; i < n ; i++ )
        {
            points.push_back( std::make_pair( i, 1 ) );
        }

    }
    else if ( edge + interior >= n )
    {
        // Do all rows with some edges
        for ( int i = 0; i < n ; i ++ )
            points.push_back( std::make_pair(i, ( i > n - 3 || i < 3 ) ? -1 : 1 ) ) ;
    }
    else
    {
        std::set<int> interior_points;
        while ( interior_points.size() < interior )
        {
            int index = rand() % n ;
            if (! ( index < edge/2.0 || index > n - edge/2  ) )
            {
                interior_points.insert(index);
            }
        }
        for ( auto it : interior_points )
        {
            points.push_back( std::make_pair( it, 1 ) ) ;
        }
        for ( int i = 0; i < edge; i++ )
        {
            points.push_back( std::make_pair( i, -1 ) );
            points.push_back( std::make_pair( n-i-1, -1 ) );
        }
    }

    //Sort the edges in order of rows.
    std::sort(points.begin(), points.end());


    // This is to remove dups and check for edge/mid combos. This is really a tiny problem problem. For
    // big matricies, this wont do anything.
    int last = 0;
    for ( int i = 1; i < points.size() ; i++ )
    {
        if ( points[i] == points[last] )
            continue; // This is a duplicate so get rid of it.
        else if ( points[i].first == points[last].first )
            points[last].second = 1; // This is a mid point and a edge point so treat it as a mid point
        else
            points[++last] = points[i]; // no match so its the first of its kind.
    }
    points.resize(last+1);
    MapDMNumberingToNaturalNumbering(points);
    return 0;

}

int
DefaultPetscTestingSpace::GetJacobianColumn(Mat J,
        std::pair< int, int > &point,
        Vec *s )
{
    Vec basis, extra; // The basis vector we will multiply by
    MatCreateVecs( J, &basis, &extra ) ; // make basis compatible with J.
    VecDuplicate(basis, s);

    PetscInt high, low;
    VecGetOwnershipRange(basis, &low, &high );

    VecSet(basis,0.0);
    VecSetValue( basis, point.first, (PetscScalar) 1.0 , INSERT_VALUES );

    VecAssemblyBegin(basis);
    VecAssemblyEnd(basis);
    MatMult( J, basis, *s );
    VecDestroy(&basis);
    VecDestroy(&extra);
    return 0;
}

int
DefaultPetscTestingSpace::GetJacobianSamplePoints(Mat J ,
        int edge,
        int interior,
        std::vector< std::pair< int, int > >&points )
{
    PetscInt n,m;
    MatGetSize(J, &m, &n ) ;
    if ( edge < 0 && interior < 0 )
        GetSamplePoints( n , edge, interior, points );
    else if ( edge < 0 || interior < 0 )
        GetSamplePoints( n, n, n, points ) ;
    else
        GetSamplePoints( n, edge, interior,points ) ;
    return 0;
}

void
DefaultPetscTestingSpace::MPI_FeatureReduce(void *invec,
        void *inoutvec,
        int *len,
        MPI_Datatype *datatype)
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

    return ;
}

int
DefaultPetscTestingSpace::ExtractJacobianFeatures(Mat J,
        int edge,
        int interior,
        std::map< std::string, double > &fnames )
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


#if NNZ
    double nnz = 0;
#endif

#if ABSOLUTENONZEROSUM
    double abs_nz_sum = 0.0; // Sum over all non zeros in the middle block
#endif

#if TRACE
    double trace = 0.0;     // trace
#endif

#if ABSOLUTETRACE
    double abs_trace = 0.0; // absolute value of the trace;
#endif

#if FROBENIUSNORM
    double fro_norm  = 0.0 ; // all middle points
#endif

#if DIAGONALAVERAGE
    double diag_ave = 0.0; // average of abs value of diagonal
#endif

#if AVERAGEDIAGONALDISTANCE
    double ave_diag_dist = 0.0; // average distance of non-zeros from the diagonal -- middle points all
    double ave_diag_dist_denom = 0 ;
#endif

#if LOWERBANDWIDTH
    double blower = 0;
#endif

#if UPPERBANDWIDTH
    double bupper = 0;
#endif

#if LOWERBANDWIDTH || UPPERBANDWIDTH
    int bwidth = 0;
#endif

#if DIAGONALMEAN
    double diag_mean = 0.0;
#endif

#if DIAGONALSIGN
    double diag_sign = 3;
#endif

#if DIAGONALNONZEROS || DIAGONALAVERAGE
    double diag_nnz = 0;
#endif

#if ROWVARIANCE
    double row_var = 0.0;
#endif

#if MAXNONZEROSPERROW
    std::vector< int > max_nnz_row(npoints); // nnz_max_rows uses sample points only from the middle only
#endif

#if MINNONZEROSPERROW
    //Set all these values to n+1 to start with. Since this is impossible ( i.e, more nz than rows ), we can use that to test weather
    //a row was accessed on this processor.
    std::vector< int > min_nnz_row(npoints, n+1); // nnz_max_rows uses sample points only from edges only
#endif

#if ONENORM
    std::vector<double> one_norm(npoints); // This is the abs sum of column for all samples this one needs mpi sum then max
#endif

#if INFINITYNORM
    std::vector<double> inf_norm(npoints); // All sample points , mid and edges
#endif

#if COLDIAGONALDOMINANCE
    std::vector<double> col_diag_dom(npoints);
#endif

#if ROWDIAGONALDOMINANCE
    // Row diagonal dominance uses a MPI over the definition for each row, with diagonal dominance only
    // true if every row is negative. So, to make sure this works in parallel, make every row negative,
    //  then set the sum to zero when this row is first accessed.
    std::vector<double> row_diag_dom(npoints,-1.0);
#endif

#if ROWVARIANCE
    std::vector<double> row_squared_sum(npoints);  // THese are for the variences s^2 _i = ( E(x^2) - n ( E x )^2 )/(n-1)
    std::vector<double> row_sum(npoints);
#endif

#if COLUMNVARIANCE
    std::vector<double> col_squared_sum(npoints);
    std::vector<double> col_sum(npoints);
#endif

#if DIAGONALSIGN
    int neg(0), pos(0),zero(0);
#endif

#if SYMMETRY
    std::vector<double> sample_data(npoints * npoints);  // all zeros, set row*npoints + column
#endif

    int num_t1 = 0;
    int num_samples = 0;
    for ( int i = 0; i < npoints; i++ )
    {
        if ( i % 100 == 0 ) std::cout << " Starting row " << i << " of " << npoints <<std::endl;

        Vec vec;
        GetJacobianColumn( J, points[i], &vec );

        VecGetArrayRead( vec, &array );
        VecGetOwnershipRange( vec, &low, &high ) ;
        type = points[i].second;
        if ( type == 1 ) nsamples_mid++;

        k = 0;
        while ( points[k].first < low && k < points.size() ) k++;   // Find the first sample point.
        for ( int j = low; j < high; j++ )
        {
            val = array[j-low];

            aval = fabs(val);


            // Is this a Sample point ?
            if ( k < points.size() && j == points[k].first )
            {
                sample=true;
                k++;
#if MINNONZEROSPERROW
                if (min_nnz_row[k-1] == n+1)
                {
                    min_nnz_row[k-1] = 0;
                }
#endif

#if ROWDIAGONALDOMINANCE
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

#if NNZ
                    nnz++;
#endif

#if ABSOLUTENONZEROSUM
                    abs_nz_sum += aval;
#endif

#if FROBENIUSNORM
                    fro_norm += aval*aval;
#endif

#if AVERAGEDIAGONALDISTANCE
                    if ( points[i].first != j )
                    {
                        ave_diag_dist += abs( points[i].first - j ) ;
                        ave_diag_dist_denom ++;
                    }
#endif

#if COLUMNVARIANCE
                    col_squared_sum[i] += aval*aval;
                    col_sum[i] += val;
#endif

                    if (sample)
                    {

#if ROWVARIANCE
                        row_squared_sum[k-1] += aval*aval;
                        row_sum[k-1] += val;
#endif

#if MAXNONZEROSPERROW
                        max_nnz_row[k-1]++;
#endif
                    }
                }
            }

#if MINNONZEROSPERROW // Count non-zeros in the sample rows. 
            if (aval > NONZEROTOLERANCE && sample ) min_nnz_row[k-1]++;
#endif

            // These are features that use both mid and edge points
#if ONENORM
            one_norm[i] += aval;
#endif

#if ( !(HAVE_DIAGONAL) )
            if (points[i].first == j )
            {
#if DIAGONALSIGN
                if ( val < -1.0*NONZEROTOLERANCE ) neg = 1;
                else if ( aval < NONZEROTOLERANCE ) zero = 1;
                else pos = 1;
#endif

                if ( sample )
                {
#if TRACE
                    trace += val;
#endif

#if ABSOLUTETRACE
                    abs_trace += aval;
#endif

#if DIAGONALNONZEROS || DIAGONALAVERAGE
                    if ( aval > NONZEROTOLERANCE ) diag_nnz++;
#endif
                }
            }
#endif

#if COLDIAGONALDOMINANCE
            col_diag_dom[i] += ( ( points[i].first == j ) ? -1 : 1 )*aval;
#endif

            if (sample)
            {
#if INFINITYNORM
                inf_norm[k-1] += aval; // points[k-1] = j, so save it in k-1.
#endif

#if ROWDIAGONALDOMINANCE
                row_diag_dom[k-1] += ( ( points[i].first == j ) ? -1 : 1 )*aval ;
#endif

#if SYMMETRY
                sample_data[(k-1)*npoints + i ] = val;
#endif
            }

#if LOWERBANDWIDTH || UPPERBANDWIDTH
            bwidth = points[i].first - j ;
            if ( aval > NONZEROTOLERANCE )
            {
#if UPPERBANDWIDTH
                if ( bwidth  > bupper ) bupper = bwidth ;
#endif
#if LOWERBANDWIDTH
                if ( -bwidth > blower ) blower = -bwidth ;
#endif
            }
#endif // LOWERBANDWIDTH || UPPERBANDWIDTH
        }
        VecRestoreArrayRead( vec , &array ) ;
        VecDestroy(&vec);
    }

#if ( DIAGONAL && HAVE_DIAGONAL )
    Vec diag, extra;
    MatCreateVecs( J, &diag, &extra ) ;
    MatGetDiagonal( J, diag );
    VecGetArrayRead( diag, &array ) ;
    VecGetOwnershipRange( diag, &low, &high );
    for ( int i = low; i < high; i++ )
    {
        val = array[i-low];
        aval = abs(val);

#if TRACE
        trace += val;
#endif

#if ABSOLUTETRACE
        abs_trace += aval;
#endif

#if DIAGONALNONZEROS || DIAGONALAVERAGE
        if ( aval > NONZEROTOLERANCE ) diag_nnz++;
#endif

#if DIAGONALSIGN
        if ( val < 0 ) neg = 1;
        else if (val == 0 ) zero = 1;
        else pos = 1;
#endif
    }
    VecRestoreArrayRead( diag, &array );
    VecDestroy(&diag);
    VecDestroy(&extra);
#endif

#if ROWVARIANCE
    row_var = 0.0;
    for ( auto i = 0; i < npoints; i++ )
        row_var += ( 1./ ( nsamples_mid*(nsamples_mid-1) ) )* ( row_squared_sum[i] - (1.0/nsamples_mid) * row_sum[i] * row_sum[i] );
#endif

    // Finish it
    int feature_count  = FEATURECOUNT(npoints);
    double *features = (double*) malloc( sizeof(double) * ( feature_count  ) );
    double *rfeatures = (double*) malloc( sizeof(double) * ( feature_count ) ) ;
    int c = 0;

#if NNZ
    features[c++] =  nnz * n / (double) nsamples_mid ;
#endif

#if AVGNONZEROSPERROW
    features[c++] = nnz / (double) nsamples_mid ;
#endif

#if ABSOLUTENONZEROSUM
    features[c++] = abs_nz_sum * n / (double) nsamples_mid;
#endif

#if TRACE
#if HAVE_DIAGONAL
    features[c++] = trace;
#else
    features[c++] = trace * n / (double) npoints  ;
#endif
#endif

#if ABSOLUTETRACE
#if HAVE_DIAGONAL
    features[c++] = abs_trace;
#else
    features[c++] = abs_trace * n / (double) npoints  ;
#endif
#endif

#if DIAGONALMEAN
#if HAVE_DIAGONAL
    features[c++] = trace / (double) n ;
#else
    features[c++] = trace / (double) npoints ;
#endif
#endif

#if (DIAGONALNONZEROS || DIAGONALAVERAGE )
#if HAVE_DIAGONAL
    features[c++] = diag_nnz ;
#else
    features[c++] = diag_nnz * n / (double) npoints ;
#endif
#endif

#if DIAGONALAVERAGE
    features[c++] = abs_trace;
#endif

#if FROBENIUSNORM
    features[c++] = fro_norm * n * n / ( nsamples_mid * nsamples_mid ) ;
#endif

#if ROWVARIANCE
    features[c++] = row_var;
#endif

#if AVERAGEDIAGONALDISTANCE
    features[c++] = ave_diag_dist ;
    features[c++] = ave_diag_dist_denom ;
#endif

    /*  Start the MPI MAX FEATURES. *****************************************/
#if INFINITYNORM
    features[c++] = *std::max_element(inf_norm.begin(),inf_norm.end());
#endif

#if MAXNONZEROSPERROW
    features[c++] = *max_element(max_nnz_row.begin(), max_nnz_row.end());
#endif

#if DIAGONALSIGN
    features[c++] = neg;
    features[c++] = pos;
    features[c++] = zero;
#endif

#if LOWERBANDWIDTH
    features[c++] = blower;
#endif

#if UPPERBANDWIDTH
    features[c++] = bupper;
#endif

#if ROWDIAGONALDOMINANCE
    features[c++] = *std::max_element(row_diag_dom.begin(), row_diag_dom.end());
#endif

    /* Start the MPI MIN Features ********************************************/
#if MINNONZEROSPERROW
    features[c++] = *min_element(min_nnz_row.begin(), min_nnz_row.end());
#endif

    /* Start the column sum features ******************************************/
#if ONENORM
    for ( auto &it : one_norm ) features[c++] = it;
#endif

#if COLDIAGONALDOMINANCE
    for ( auto &it : col_diag_dom ) features[c++] = it;
#endif

#if COLUMNVARIANCE
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
#if NNZ
        fnames.insert( std::make_pair( "nnz",                  rfeatures[c++] ) );
#endif

#if AVGNONZEROSPERROW
        fnames.insert( std::make_pair( "AvgNonzerosPerRow",    rfeatures[c++]  ) );
#endif

#if ABSOLUTENONZEROSUM
        fnames.insert( std::make_pair( "AbsoluteNonZeroSum",   rfeatures[c++]  ) );
#endif

#if TRACE
        fnames.insert( std::make_pair( "Trace",                rfeatures[c++]  ) );
#endif

#if ABSOLUTETRACE
        fnames.insert( std::make_pair( "AbsoluteTrace",        rfeatures[c++]  ) );
#endif

#if DIAGONALMEAN
        fnames.insert( std::make_pair( "DiagonalMean",         rfeatures[c++]  ) );
#endif

#if (DIAGONALAVERAGE || DIAGONALNONZEROS  )
        fnames.insert( std::make_pair( "DiagonalNonZeros",     rfeatures[c++]  ) );
#endif

#if DIAGONALAVERAGE
        fnames.insert( std::make_pair( "DiagonalAverage",      rfeatures[c++]/fnames["DiagonalNonZeros"] ) );
#endif

#if FROBENIUSNORM
        fnames.insert( std::make_pair( "FrobeniusNorn",        sqrt(rfeatures[c++])) );
#endif

#if ROWVARIANCE
        fnames.insert( std::make_pair( "RowVariance",          rfeatures[c++] ) );
#endif

#if AVERAGEDIAGONALDISTANCE
        double sum = rfeatures[c++];
        double denom = rfeatures[c++];
        fnames.insert( std::make_pair( "AvgDiagonalDistance" , sum/denom ) );
#endif

#if INFINITYNORM
        fnames.insert( std::make_pair( "InfinityNorm",         rfeatures[c++] ) );
#endif

#if MAXNONZEROSPERROW
        fnames.insert( std::make_pair( "MaxNonZeroPerRow",     rfeatures[c++] ) );
#endif

#if DIAGONALSIGN
        // Perform all the buffer post-processing ..
        neg = rfeatures[c++];
        pos = rfeatures[c++];
        zero = rfeatures[c++];
        if ( neg && pos ) diag_sign = -3;
        else if ( neg && zero ) diag_sign = -1;
        else if ( pos && zero ) diag_sign = -1;
        else if ( neg ) diag_sign = -2;
        else if ( pos ) diag_sign = 2;
        else if ( zero) diag_sign = 0;
        fnames.insert( std::make_pair( "DiagonalSign",         diag_sign ) );
#endif

#if LOWERBANDWIDTH
        fnames.insert( std::make_pair( "lowerBandwidth",       rfeatures[c++] ) );
#endif

#if UPPERBANDWIDTH
        fnames.insert( std::make_pair( "upperBandwidth",       rfeatures[c++] ) );
#endif

#if ROWDIAGONALDOMINANCE
        fnames.insert( std::make_pair( "RowDiagonalDominance", ( rfeatures[c++] < 0 ) ? 1 : 0 ) );
#endif

#if MINNONZEROSPERROW
        fnames.insert( std::make_pair( "MinNonZeroPerRow",     rfeatures[c++] ) );
#endif

#if ONENORM
        double fone_norm = -1;
        for ( int i = 0; i < npoints; i++ )
        {
            t = rfeatures[c++];
            if ( t > fone_norm ) fone_norm = t;
        }
        fnames.insert( std::make_pair( "OneNorm", fone_norm ) );
#endif

#if COLDIAGONALDOMINANCE
        double ct;
        double fccd = 1;
        for ( int i = 0; i < npoints; i++ )
        {
            t = rfeatures[c++];
            if ( t > 0 ) fccd = 0;
        }
        fnames.insert( std::make_pair( "ColumnDiagonalDominance", fccd ) );
#endif

#if COLUMNVARIANCE
        double cvar = 0;
        for ( int i = 0; i < npoints; i++ )
        {
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

#if SYMMETRICITY
            sym = ( sym && fabs( value - svalue ) < 1e-13 );
#endif

#if NONZEROPATTERNSYMMETRY
            nzsym = ( nzsym && (  ( fabs(value) > 1e-13 && fabs(svalue) > 1e-13 ) || ( fabs(value) < 1e-13 && fabs(svalue) < 1e-13 ) ) );
#endif

#if SYMMETRICINFINITYNORM
            sinfnorm[row] += fabs( 0.5 * ( svalue + value ) );
#endif

#if ANTIYSYMMETRICINFINITYNORM
            asinfnorm[row] += fabs( 0.5 * ( svalue - value ) );
#endif

#if (SYMMETRICFROBENIUSNORM || ANTISYMMETRICFROBENIUSNORM)
            if ( points[col].second == 1 )
            {
#if SYMMETRICFROBENIUSNORM
                sfronorm += ( 0.5 * ( svalue + value ) ) * ( 0.5 * ( svalue + value ) ) ;
#endif

#if ANTISYMMETRICFROBENIUSNORM
                asfronorm += ( 0.5 * ( svalue - value ) ) * ( 0.5 * ( svalue - value ) ) ;
#endif
            }
#endif

        }
#if SYMMETRICITY
        fnames.insert( std::make_pair( "Symmetricity", sym  ) );
#endif

#if NONZEROPATTERNSYMMETRY
        fnames.insert( std::make_pair( "NonZeroPatternSymmetry", nzsym ) );
#endif

#if SYMMETRICINFINITYNORM
        fnames.insert( std::make_pair( "SymmeticInfinityNorm", *std::max_element(sinfnorm.begin(), sinfnorm.end()) ));
#endif

#if ANTISYMMETRICINFINITYNORM
        asfronorm = n * sqrt(asfronorm) / (double) nsamples_mid ;
        fnames.insert( std::make_pair( "AntiSymmetricFrobeniusNorm" , asfronorm ) );
#endif

#if SYMMETRICFROBENIUSNORM
        sfronorm = n * sqrt(sfronorm) / (double) nsamples_mid ;
        fnames.insert( std::make_pair( "SymmetricFrobeniusNorm", sfronorm   ) );
#endif

#if ANTISYMMETRICFROBENIUSNORM
        asfronorm = n * sqrt(asfronorm) / (double) nsamples_mid ;
        fnames.insert( std::make_pair( "AntiSymmerticInfinityNorm", *std::max_element(asinfnorm.begin(),asinfnorm.end())) );
#endif

#endif  // SYMMETRY  
    }

    free(features);
    free(rfeatures);
    return 0;
}

}
#endif 
