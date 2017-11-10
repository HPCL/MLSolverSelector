
/* \file MLWaff.C
 * \brief Implimentation of the ML base class using Waffles ML package. 
 **/

#ifdef WITH_WAFFLES

#include "MLWaff.h"


using namespace MachineLearning; 

_ML_Waffles::_ML_Waffles() : _SS_MachineLearning() 
{                            
 features = NULL;
 labels = NULL;
 model = NULL;
 fmodel = NULL;

 Add("method","KNN","ML method to use"); /* the ML algorithm to use ( KNN by defalt ) */
}

_ML_Waffles::~_ML_Waffles()
{ 
 if ( model ) delete model;
 if ( fmodel ) delete fmodel; 
 if ( features) delete features;
 if ( labels ) delete labels;
}

_SS_ErrorFlag
_ML_Waffles::ImportData( std::shared_ptr<GClasses::GMatrix> &matrix,
                         std::vector< int > &labelsdim,  
                         std::vector< std::pair< std::string, std::string > > &sset )
{ 
  
 std::vector< std::vector < std::string >> data;
 std::vector < std::string > column_names;
 std::vector< int > feature_or_label; 
 database->ImportData( sset  , data , column_names, feature_or_label);
 
 GClasses::GArffRelation *pRelation = new GClasses::GArffRelation();
 pRelation->setName("SolverSelecter");

 for (unsigned int i = 0; i < column_names.size(); i++) 
 { 
    if ( feature_or_label[i] == 0 )
    {
        pRelation->addAttribute( column_names[i].c_str(), 0, NULL );
    }
    else
    {
       labelsdim.push_back(i);
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
      (*matrix)[row_count][0] = std::stoi( row[0] );
      for ( int i = 1; i < size ; i++  )
      {
	        size_t vals = pRelation->valueCount(i);
          if ( vals == 0 )
            (*matrix)[row_count][i] = std::atof(row[i].c_str());
          else 
            (*matrix)[row_count][i] = (double) pRelation->findEnumeratedValue(i,row[i].c_str()); 		      
      }
      row_count++;
	}
  return _SS_error_flag ;
}


_SS_ErrorFlag 
_ML_Waffles::TrainSystem( std::vector< std::pair< std::string, std::string > > &sset )
{

  
  std::shared_ptr< GClasses::GMatrix > matrix(nullptr) ;
 std::vector< int > labelsdim;
 ImportData( matrix, labelsdim, sset ) ; 

 const GClasses::GRelation &x = matrix->relation();
 for ( unsigned int i = 0; i < x.size(); i++ )
 {
    attrs.push_back( x.attrNameStr(i) ) ;
 }
 
 /* Split by features and labels */
 GClasses::GDataColSplitter splitter(*matrix, (std::size_t) labelsdim.size() );
 if ( features ) { delete features; } features = new GClasses::GMatrix();
 if ( labels ) { delete labels ; } labels   = new GClasses::GMatrix(); 
 
 *features = splitter.features();
 *labels   = splitter.labels();
  
 features_c = features->cols();
 features_r = features->rows();
 labels_c = labels->cols();
 labels_r = labels->rows();
  
 /* Build the model -- Models all seem to inheret from the supervised learner */
 std::string method;
 Get("method",method);
 if ( method == "DescisionTree" )
    model = new GClasses::GDecisionTree();
 else if ( method == "RandomForest" ) 
    model = new GClasses::GRandomForest(50,2);
 else if ( method == "KNN")
 {
     model = new GClasses::GKNN();
     GClasses::GKNN *m = ( GClasses::GKNN*) model;
     m->autoTune( *features, *labels );   
 }
 else if ( method == "Basyian" )
    model = new GClasses::GNaiveBayes();
 else
    model = new GClasses::GBaselineLearner();

 /* Filter the results for any unsupported data types */
 fmodel = new GClasses::GAutoFilter( model, false );     
 
 fmodel->train( *features , *labels );
 return _SS_error_flag;

}

_SS_ErrorFlag 
_ML_Waffles::Predict( _SS_Features &afeatures, const int &solver_hash , std::vector<bool> &good  )
{
  GClasses::GVec pattern(features_c);
  const GClasses::GRelation &x = features->relation();
 
  
  std::map< std::string, double > feat_map;
  afeatures.Get(feat_map);
  int i = 0;
  while ( i < x.size() && x.attrNameStr(i) != "HASH" ) i++;
  pattern[i] = solver_hash ;
  
  for ( auto it : feat_map )
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

#endif

