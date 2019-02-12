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

//change from loop to strng and double
ErrorFlag
DatabaseInterface::AddClassification( int  &urow ,
                                      std::map< std::string, double > &cvalues,
                                      std::string CName ,UpdateStatus status )
{
    std::map<std::string, double> mvalues;
    std::map<std::string, bool >  bvalues;
    GetUniqueMeasurement( urow , mvalues );

    for ( auto &it : cvalues )
    {
        auto m = mvalues.find(it.first);
        if ( m == mvalues.end() || m->second < 0.0 )
        {
            bvalues[CName] = false;
        }
        else 
        {
            bvalues[CName] = ( m->second < it.second ) ? true : false ;
        }
    }
    //printf("urow %d\n",urow);
    AddClassificationToDatabase( urow, bvalues, UpdateStatus::REPLACE  );
    return 1;
}

ErrorFlag
DatabaseInterface::ClassifySolvers( std::string MName, double MValue, std::string CName  )
{

    std::map< int, std::vector< int > > matrixmap;

    GetMatrixToUniqueMap( matrixmap, -1 ) ;

    for ( auto matrix : matrixmap )
    {
        std::map<std::string, double> cvalues;
        for ( auto solver : matrix.second )
        {
            std::map<std::string, double> mvalues;
            GetUniqueMeasurement( solver , mvalues ) ;

                auto m = mvalues.find(MName);
                auto c = cvalues.find(MName);
                
                    
                if ( m != mvalues.end() )
                {
                    
                    double cutoff = ( 1.0 + MValue )*m->second ;
                    if ( c == cvalues.end() )
                        cvalues[MName] = ( m->second < 0 ) ? -1 : cutoff;
                    else if ( cutoff < c->second && m->second > 0 )
                    {
                        c->second = cutoff ;
                    }
                }
            
        }


        
        for ( auto solver : matrix.second )
        {
         //   for ( auto &it : cvalues )
                AddClassification( solver, cvalues, CName,UpdateStatus::REPLACE ) ;
        }
    }
    return 1;
}

ErrorFlag
DatabaseInterface::GetMachineLearningData(std::vector <std::string>&CNames,std::vector< int > &row_ids,
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
    
    if ( CNames.size() > 0 ) {
    
        std::vector<std::string> classification_temp;
        GetClassificationLabels(classification_temp);
   
        for ( auto want : CNames ) {
            auto f = std::find(classification_temp.begin(), classification_temp.end(), want); //find(want);
            if  ( f != classification_temp.end() )
            {
                classification_labels.push_back(want);
            }
            else {
                printf("CName %s not in database" , want.c_str());
            }
        }
        if ( classification_labels.size() == 0 ) {
            printf(" No classifications Found \n");
            std::abort();
        }
    } else {
         GetClassificationLabels(classification_labels);
    }

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


