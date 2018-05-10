#ifndef _ML_WAFFLES_H
#define _ML_WAFFLES_H

/** \file MLWaff.h
 * \brief Implimentations for the ML waffles class. */
#include "MachineLearning.h"
#include "GClasses/GMatrix.h"
#include "GClasses/GLearner.h"
#include "GClasses/GKNN.h"
#include "GClasses/GDecisionTree.h"
#include "GClasses/GNaiveBayes.h"
#include "GClasses/GTransform.h"

/**
 * The _ML_Waffles class does the actual ML for the Machine learning class. Need
 * to change the name to _SS_Waffles for consistency
 **/
class _ML_Waffles : public _SS_MachineLearning
{
public:

    std::string method; /**< name of the method */
    GClasses::GSupervisedLearner *model;    /**< the machine learning model */
    GClasses::GAutoFilter        *fmodel;   /**< a filtered machine learning model */
    GClasses::GMatrix            *features;  /**< waffles matrix containing the features */
    GClasses::GMatrix            *labels ;  /**< waffles matrix containing the labels */


    int features_c; /**< number of feature columns */
    int features_r; /**< number of feature rows */
    int labels_r;   /**< number of label rows */
    int labels_c;   /**< number of label columns */

    /** Default constructer */
    _ML_Waffles() {
        features = NULL;
        labels = NULL;
        model = NULL;
        fmodel = NULL;
        method = "KNN";

    }


    /** Default destructor */
    ~_ML_Waffles() {
        if ( model ) delete model;
        if ( fmodel ) delete fmodel;
        if ( features) delete features;
        if ( labels ) delete labels;
    }


    /**
     * Train the chosen model. The training set can be restricted to only certain solver/preconditioner
     * pairs. If sset.size() > 0, only this pairs listed in sset are used for training, and hence, for
     * solving the problem. Future versions might allow for restriction based on parameters as well.
     **/
    _SS_ErrorFlag TrainSystem(  ) override {
        
        std::vector< std::string > fnames ;
        database->GetFeatureNames(fnames);
        BuildModel( method, fnames );
        fmodel->train( *features , *labels );
        return _SS_error_flag;
    }
    
    /**
     * Classify a matrix based on its _SS_Features and return chosen solver
     **/
    _SS_ErrorFlag Classify( _SS_features_map &features /**< the feature set of the matrix */,
                            _SS_Solver &solver /**< output, a (hopefully) "good" solver for the problem */) override
    {

        /* Classify the system, choose the best solver, get the parameter string for
         * the solver, and initialize the solver with it */
        std::vector< bool > good;

        bool found_one = false;
        std::vector< int > hash_list;
        database->GetSolverHashList(hash_list);

        for ( auto hash : hash_list )
        {
            Predict( features, hash, good );
            auto it = good.begin();
            while ( it != good.end() )
            {
                if (*it) it++;
                else break;
            }
            if ( it == good.end() )
            {
                found_one = true;
                database->GetSolver( hash, solver );

                break;
            }
        }
        if ( !found_one )
        {
            std::cout<<"No Good Solvers -- Using the default \n";
        }

        return _SS_error_flag;

    }


    /**
     * Predict if a feature set is going to be a good solver.
     **/
    _SS_ErrorFlag Predict( _SS_features_map &afeatures, /**< feature set to test */
                           const int &solver_hash,    /**< hash of the solver to test */
                           std::vector<bool> &good /**< bool staing if solver is good in each category */
                         ) 
    {
        GClasses::GVec pattern(features_c);
        const GClasses::GRelation &x = features->relation();


        int i = 0;
        while ( i < x.size() && x.attrNameStr(i) != "HASH" ) i++;
        pattern[i] = solver_hash ;

        for ( auto it : afeatures )
        {
            i = 0;
            while ( i < x.size() && it.first != x.attrNameStr(i) )
            {
                i++;
            }

            if ( i == x.size() )
                std::cout << it.first << " not found \n ";
            else
            {
                pattern[i] = it.second;

            }
        }


        /* The prediction is a vector equal to the length of the labels. If we are just
         * measureing time, this is a vector fo length one */
        GClasses::GVec prediction(labels_c);

        /* Get the prediction */
        fmodel->predict( pattern ,prediction );

        good.clear();
        for ( int i = 0; i < labels_c; i++ )
        {
            good.push_back( (bool) prediction[i] );
        }

        return _SS_error_flag;
    }


    /** Extra function that converts the data into a GMatrix */
    _SS_ErrorFlag ImportData(
        std::vector < std::string > &feature_list,
        std::shared_ptr<GClasses::GMatrix> &matrix,
        int &num_labels,
        std::vector< std::pair< std::string, std::string > > &sset ) {

        std::vector< std::vector < std::string >> data;
        std::vector < std::string > column_names;
        std::vector< int > feature_or_label;
        database->ImportData( sset  , data , column_names, feature_or_label);

        GClasses::GArffRelation *pRelation = new GClasses::GArffRelation();
        pRelation->setName("SolverSelecter");

        std::vector< int > add_features, add_labels;
        for (unsigned int i = 0; i < column_names.size(); i++)
        {
            if ( feature_or_label[i] == 0 )
            {
                if ( i == 0 || std::find(feature_list.begin(),feature_list.end(), column_names[i] ) != feature_list.end() )
                {
                    pRelation->addAttribute( column_names[i].c_str(), 0, NULL );
                    add_features.push_back(i);
                }
            }
        }
        for (unsigned int i = 0; i < column_names.size(); i++)
        {
            if ( feature_or_label[i] != 0 )
            {
                add_labels.push_back(i);
                num_labels++;
                std::vector< const char *> pValues = {"0","1"};
                pRelation->addAttribute( column_names[i].c_str(), 2, &pValues);
            }
        }
        /* I believe GMatrix takes control of pRelation, and deletes it */
        matrix = std::make_shared<GClasses::GMatrix>(pRelation);

        matrix->newRows( data.size() );
        int row_count = 0;
        for ( auto row : data )
        {
            int size = matrix->cols();
            int col_count = 0;
            for ( auto it : add_features )
            {
                if ( it == 0 )
                    (*matrix)[row_count][col_count] = std::stoi( row[it].c_str() );
                else
                    (*matrix)[row_count][col_count] = std::atof(row[it].c_str());
                col_count++;
            }
            for ( auto it : add_labels )
            {
                (*matrix)[row_count][col_count] = (double) pRelation->findEnumeratedValue(col_count,row[it].c_str());
                col_count++;
            }
            row_count++;
        }

        return _SS_error_flag ;
    }



    _SS_ErrorFlag CrossValidate( std::string algorithm, std::vector< std::string > &features_list ) override {


        std::cout << " ----------   Performing Cross Validation -------------- \n";
        std::cout << " --- Features -- " ;
        for ( auto it : features_list ) std::cout << it << " -- ";
        std::cout << "\n\n";

        // Parse options
        unsigned int seed = (unsigned int)time(NULL);
        int reps = 5;
        int folds = 5;
        bool succinct = true;

        if(reps < 1)
            throw GClasses::Ex("There must be at least 1 rep.");
        if(folds < 2)
            throw GClasses::Ex("There must be at least 2 folds.");
        algorithm = "KNN";
        BuildModel( algorithm , features_list );
        fmodel->rand().setSeed(seed);

        // Test
        std::cout.precision(8);
        double sae;
        double sse = fmodel->repValidate(*features, *labels, reps, folds, &sae, succinct ? NULL : CrossValidateCallback, fmodel);
        if(succinct)
        {
            std::cout << "Misclassification : " << GClasses::to_str(sse / features->rows()) << std::endl ;
            std::cout << "Predictive Acc : " << GClasses::to_str(1.0 - (sse / features->rows())) << std::endl;;
        }
        else
        {
            if(labels->cols() == 1 && labels->relation().valueCount(0) > 0)
            {
                std::cout << "Misclassification rate: " << GClasses::to_str(sse / features->rows()) << "\n";
                std::cout << "Predictive accuracy: " << GClasses::to_str(1.0 - (sse / features->rows())) << "\n";
            }
            else
            {
                std::cout << "Sum absolute error: " << GClasses::to_str(sae) << "\n";
                std::cout << "Mean absolute error: " << GClasses::to_str(sae / features->rows()) << "\n";
                std::cout << "Sum squared error: " << GClasses::to_str(sse) << "\n";
                std::cout << "Mean squared error: " << GClasses::to_str(sse / features->rows()) << "\n";
                std::cout << "Root mean squared error: " << GClasses::to_str(sqrt(sse / features->rows())) << "\n";
            }
        }

        std::cout << " ------------ Finished Cross Validation with " << algorithm << " ------------------- \n\n" ;

    }




    _SS_ErrorFlag BuildModel(  std::string method, std::vector< std::string > &features_list) {
        /* This is a place holder to support restricting the solver set. */
        std::vector< std::pair< std::string, std::string > > sset;

        std::shared_ptr< GClasses::GMatrix > matrix(nullptr) ;
        int labelsdim;
        ImportData( features_list, matrix, labelsdim, sset ) ;

        const GClasses::GRelation &x = matrix->relation();

        /* Split by features and labels */
        GClasses::GDataColSplitter splitter(*matrix, (std::size_t) labelsdim );
        if ( features ) {
            delete features;
        }
        features = new GClasses::GMatrix();
        if ( labels ) {
            delete labels ;
        }
        labels   = new GClasses::GMatrix();

        *features = splitter.features();
        *labels   = splitter.labels();


        features_c = features->cols();
        features_r = features->rows();
        labels_c = labels->cols();
        labels_r = labels->rows();

        /* Build the model -- Models all seem to inheret from the supervised learner */
        if ( method == "DescisionTree" )
            model = new GClasses::GDecisionTree();
        else if ( method == "RandomForest" )
            model = new GClasses::GRandomForest(50,2);
        else if ( method == "KNN")
        {
            model = new GClasses::GKNN();
            GClasses::GKNN *m = ( GClasses::GKNN*) model;
        }
        else if ( method == "Basyian" )
            model = new GClasses::GNaiveBayes();
        else
            model = new GClasses::GBaselineLearner();

        /* Filter the results for any unsupported data types */
        if ( fmodel ) {
            delete fmodel ;
        }
        fmodel = new GClasses::GAutoFilter( model, false );

        return _SS_error_flag;
    }


    static void CrossValidateCallback(void* pSupLearner, size_t nRep, size_t nFold, double foldSSE, size_t rows) {
        std::cout << "Rep: " << nRep << ", Fold: " << nFold <<", Mean squared error: " << GClasses::to_str(foldSSE / rows) << "\n";
    }

};

#endif

