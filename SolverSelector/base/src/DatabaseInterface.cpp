#include "DatabaseInterface.h"

namespace SolverSelecter
{

DatabaseInterface::DatabaseInterface() : SSBase("DATABASE")
{
  initialized = false;
}

ErrorFlag
DatabaseInterface::Finalize()
{
    return 0;
}

ErrorFlag
DatabaseInterface::AddRow(Solver &solver,
                          std::map< std::string, double > &features,
                          std::map<std::string, double> &measurements)
{
    //int solver_hash = AddSolverToDatabase( solver, status  ) ;
    //int matrix_hash = AddMatrixToDatabase( features, status );
    //row_hash = AddRowToDatabase( row_hash, solver_hash, matrix_hash, measurements, status);
    
    int row_hash = AddToDatabase(solver, features, measurements, status );

    
    return row_hash;
}

ErrorFlag
DatabaseInterface::AddClassification( int  &urow ,
                                      std::map< std::string, double > &cvalues,
                                      UpdateStatus status )
{
    std::map<std::string, double> mvalues;
    std::map<std::string, bool >  bvalues;
    GetUniqueMeasurement( urow , mvalues );

    for ( auto &it : cvalues )
    {
        auto m = mvalues.find(it.first);
        if ( m == mvalues.end() || m->second < 0.0 )
        {
            bvalues[it.first] = false;
        }
        else
        {
            bvalues[it.first] = ( m->second < it.second ) ? true : false ;
        }
    }
    AddClassificationToDatabase( urow, bvalues, UpdateStatus::REPLACE  );
    return 1;
}

ErrorFlag
DatabaseInterface::ClassifySolvers(std::map< std::string ,
                                   double > &bvalues,
                                   int just_matrix_in_row )
{

    std::map< int, std::vector< int > > matrixmap;

    GetMatrixToUniqueMap( matrixmap, just_matrix_in_row ) ;

    for ( auto matrix : matrixmap )
    {
        std::map<std::string, double> cvalues;
        for ( auto solver : matrix.second )
        {
            std::map<std::string, double> mvalues;
            GetUniqueMeasurement( solver , mvalues ) ;

            for ( auto meas : bvalues )
            {
                auto m = mvalues.find(meas.first);
                auto c = cvalues.find(meas.first);
                if ( m != mvalues.end() && m->second > 0 )
                {
                    double cutoff = ( 1.0 + meas.second )*m->second ;
                    if ( c == cvalues.end() )
                        cvalues[meas.first] = cutoff;
                    else if ( cutoff < c->second )
                    {
                        c->second = cutoff;
                    }
                }
            }
        }
        for ( auto solver : matrix.second )
        {
            for ( auto &it : cvalues )
                AddClassification( solver, cvalues, UpdateStatus::REPLACE ) ;
        }
    }
    return 1;
}

ErrorFlag
DatabaseInterface::GetMachineLearningData(std::vector< int > &row_ids,
        std::vector< int > &solvers_labels,
        std::vector< std::string > &features_labels ,
        std::vector< std::string > &classification_labels,
        std::vector< int > &solvers_data,
        std::vector< std::vector< double > > &feature_data,
        std::vector< std::vector< bool > > &classification_data )
{


    std::vector< int > rows;
    GetUniqueHashList( rows );  // List of every row in the database.

    GetFeatureLabels(features_labels);
    GetClassificationLabels(classification_labels);
    GetUniqueSolverList(solvers_labels);

    for ( auto it : rows )
    {

        row_ids.push_back( it );

        std::map< std::string, double > fvalues ;
        std::map<std::string, bool > cvalues;
        std::vector< double > features;
        std::vector< bool > classi;
        int solver;

        GetUniqueSolverInRow(it, solver );
        GetUniqueFeatures( it, fvalues );
        GetUniqueClassification( it, cvalues );

        solvers_data.push_back(solver);

        for ( auto itt : features_labels )
        {
            auto f = fvalues.find(itt);
            if ( f != fvalues.end() )
                features.push_back( f->second ) ;
            else
                features.push_back( 1e300 );
        }
        feature_data.push_back(features);

        for ( auto itt : classification_labels )
        {
            auto f = cvalues.find(itt);
            if ( f != cvalues.end() )
                classi.push_back( f->second ) ;
            else
                classi.push_back( false );
        }
        classification_data.push_back(classi);
    }
    return 0;
}

}


