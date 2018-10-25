#if WITH_WAFFLES 

#include "WafflesInterface.h"

namespace SolverSelecter
{


WafflesInterface::WafflesInterface()
{
    features = NULL;
    labels = NULL;
    model = NULL;

    AddParameter("algorithm","RandomForest","Choose a machine learning model");
    AddParameter("serialized","","Filename for serialized model (empty for build at runtime");
}

WafflesInterface::~WafflesInterface()
{
    if ( model ) { delete model; model = NULL; }
    if ( features) delete features;
    if ( labels ) delete labels;
}

ErrorFlag WafflesInterface::BuildImpl()
{
    if (model == NULL ) {
      std::string serial = GetParameter("serialized"); 
      if ( serial.empty() ) 
      {
        BuildModelAtRuntime( ) ;
        database->GetUniqueSolverList(solver_hash_list);   
      }
      else
        BuildModelFromSerial(serial);
    }
    return error_flag;
}

ErrorFlag WafflesInterface::BuildModelAtRuntime()
{
    if ( model == NULL ) 
    {
      std::vector< std::string > fnames;
      database->GetFeatureLabels(fnames);
      std::string algorithm = GetParameter("algorithm");
      std::shared_ptr< GClasses::GMatrix > matrix(nullptr) ;
      int labelsdim;
      ImportData( fnames, matrix, labelsdim ) ;

      const GClasses::GRelation &x = matrix->relation();

      /* Split by features and labels */
      GClasses::GDataColSplitter splitter(*matrix, (std::size_t) labelsdim );
      if ( features )
      {
          delete features;
      }
      features = new GClasses::GMatrix();
      if ( labels )
      {
          delete labels ;
      }
      labels   = new GClasses::GMatrix();

      *features = splitter.features();
      *labels   = splitter.labels();
   
      std::string method = algorithm; //GetParameter("algorithm"); 
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

      model->train( *features , *labels );
      
      const GClasses::GRelation &xx = features->relation();
      for ( int i = 0; i < xx.size(); i++ ) features_order.push_back( xx.attrNameStr(i)) ;
      const GClasses::GRelation &xy = labels->relation();
      for ( int i = 0; i < xy.size(); i++ ) labels_order.push_back( xy.attrNameStr(i)) ;
    }
    std::cout << "FINISHED MODEL BUILDING " << std::endl;
    return error_flag;
}

ErrorFlag WafflesInterface::SerializeImpl(std::string outputfile)
{
    BuildModelAtRuntime( ) ;

    GClasses::GDom gdom;
    GClasses::GDomNode *gnode = model->serialize(&gdom);
    gnode->addField(&gdom, "method", gdom.newString(GetParameter("algorithm").c_str()));
    
    GClasses::GDomNode *ggnode = gdom.newObj(); 
    
    Solver solver;
    std::string solverstring;
    GClasses::GDomNode *slist = gdom.newList();
    GClasses::GDomNode *flist = gdom.newList();
    GClasses::GDomNode *llist = gdom.newList();
     
    for ( auto it: solver_hash_list ) 
    {
        database->GetUniqueSolver( it, solver );
        solver.GetSolverString(solverstring);
        ggnode->addField(&gdom, std::to_string(it).c_str(), gdom.newString(solverstring.c_str()) ); 
        slist->addItem(&gdom, gdom.newInt( it ));
    }
    for ( auto it : features_order ) 
      flist->addItem(&gdom, gdom.newString(it.c_str()) );
    for ( auto it : labels_order ) 
      llist->addItem(&gdom, gdom.newString(it.c_str()) );

    ggnode->addField(&gdom, "slist", slist );
    ggnode->addField(&gdom, "flist", flist );
    ggnode->addField(&gdom, "llist", llist );
    gnode->addField(&gdom, "solvers", ggnode );

    gdom.setRoot(gnode);
    gdom.saveJson(outputfile.c_str());
    return error_flag;
}

ErrorFlag WafflesInterface::BuildModelFromSerial(std::string serialized)
{
    
   if ( model == NULL ) { 
      GClasses::GDom gdom;
      gdom.loadJson(serialized.c_str());
      const GClasses::GDomNode *gnode = gdom.root();
      std::string meth = gnode->field("method")->asString();
      std::string method = GetParameter("algorithm");
      
      if ( method == "DescisionTree" )
        model = new GClasses::GDecisionTree(gnode);
      else if ( method == "RandomForest" ) {
        GClasses::GLearnerLoader ll; 
        model = new GClasses::GRandomForest(gnode, ll);
      }
      else if ( method == "KNN")
          model = new GClasses::GKNN(gnode);
      else if ( method == "Basyian" )
          model = new GClasses::GNaiveBayes(gnode);
      else
          model = new GClasses::GBaselineLearner(gnode);
      
      const GClasses::GDomNode *gg = gnode->field("solvers");
      for(GClasses::GDomListIterator li(gg->field("slist")); li.current(); li.advance())
      {
          int hash = li.current()->asInt() ; 
          solver_hash_list.push_back( hash );
          solver_hash_map[hash] = gg->field(std::to_string(hash).c_str())->asString(); 
      }
      for(GClasses::GDomListIterator li(gg->field("flist")); li.current(); li.advance())
      {
          std::string hash = li.current()->asString() ; 
          features_order.push_back( hash );
      }
      for(GClasses::GDomListIterator li(gg->field("llist")); li.current(); li.advance())
      {
          std::string hash = li.current()->asString() ; 
          labels_order.push_back( hash );
      }
  }
  return error_flag;
} 
  

ErrorFlag WafflesInterface::ClassifyImpl( features_map &afeatures /**< the feature set of the matrix */,
                                          Solver &solver /**< output, a (hopefully) "good" solver for the problem */)
{
    std::vector< bool > good;
    
    bool found_one = false;
    
    std::string serial = GetParameter("serialized"); 
    
    GClasses::GVec prediction(labels_order.size());
    GClasses::GVec pattern(features_order.size()); 

    int i = 1;
    for ( auto it : afeatures )
    {
        auto found = std::find(features_order.begin(), features_order.end(), it.first );
        if ( found != features_order.end() )
        {
          printf("\t\t %s found in model \n", it.first.c_str());
            pattern[i++] = it.second;  
        }
         else
            printf("\t\t %s not found in model \n", it.first.c_str());
    }

    for ( auto hash : solver_hash_list )
    {
        pattern[0] = hash; 
        model->predict( pattern ,prediction );
        
        bool bad = false;
        for ( int j = 0; j < prediction.size(); j++ ) 
        {
            if ( ! ((bool) prediction[j]) ) 
            {
               bad = true; 
               
            }
        }
        
        Solver tempSolver;
        database->GetUniqueSolver( hash, tempSolver );
        std::string solverT;
        tempSolver.GetSolverString(solverT);
        printf( " Solver %s  was %s \n " , solverT.c_str(), ( bad )? "bad" : "good " );

        if ( ! bad )
        {
            if (serial.empty()) {
              database->GetUniqueSolver( hash, solver );
            } else {
              solver.ParseSolverString( solver_hash_map[hash] );  
            }
            if ( ! isBaned(solver) ) {
              return 0;
            }
        }
    }
    if ( !found_one )
    {
        std::cout<<"No Good Solvers -- Using the default \n";
        solver.Clear();
    }
    return error_flag;
}

ErrorFlag WafflesInterface::ImportData(
    std::vector < std::string > &feature_list,
    std::shared_ptr<GClasses::GMatrix> &matrix,
    int &num_labels )
{

    std::vector<int> row_ids;
    std::vector<std::string> feature_labels, classification_labels;
    std::vector<std::vector< double >> feature_data;
    std::vector<std::vector< bool >> classification_data;
    std::vector<int> solvers_labels, solvers_data;
    
    GetMachineLearningData(row_ids,
                           solvers_labels,
                           feature_labels ,
                           classification_labels,
                           solvers_data,
                           feature_data,
                           classification_data );
      

    GClasses::GArffRelation *pRelation = new GClasses::GArffRelation();
    pRelation->setName("SolverSelecter");


    num_labels = classification_labels.size();

    int n = 0;
    std::vector< int > add_features, add_classis ; // features that are actually added.



    //pRelation->addAttribute("HASH", svalues.size(), &ssvalues);
    pRelation->addAttribute("HASH", 0, NULL);

    for ( auto it : feature_labels )
    {
        if ( std::find(feature_list.begin(),feature_list.end(), it.c_str() ) != feature_list.end() )
        {
            pRelation->addAttribute( it.c_str(), 0, NULL );
            add_features.push_back(n);
        }
        n++;

    }

    n = 0;
    std::vector< const char *> pValues = {"0","1"};
    for ( auto it : classification_labels )
    {
        pRelation->addAttribute( it.c_str(), 2, &pValues);
        add_classis.push_back(n);
        n++;
    }

    /* I believe GMatrix takes control of pRelation, and deletes it */
    matrix = std::make_shared<GClasses::GMatrix>(pRelation);

    matrix->newRows( row_ids.size() );
    int row_count = 0;
    int size = matrix->cols();
    for (auto it : row_ids   )
    {
        (*matrix)[row_count][0] = solvers_data[row_count] ;
        n=1;
        for ( auto it : add_features )
            (*matrix)[row_count][n++] = feature_data[row_count][it];
        for ( auto it : add_classis )
        {
            (*matrix)[row_count][n] = (double) pRelation->findEnumeratedValue(n,std::to_string(classification_data[row_count][it]).c_str());
            n++;
        }
        row_count++;
    }


    return error_flag ;
}

ErrorFlag WafflesInterface::CrossValidateImpl(int folds )
{


    std::cout << " ----------   Performing Cross Validation --------------" << GetParameter("algorithm") << " \n";

    // Parse options
    unsigned int seed = (unsigned int)time(NULL);
    int reps = 5;
    bool succinct = true;

    if(reps < 1)
        throw GClasses::Ex("There must be at least 1 rep.");
    if(folds < 2)
        throw GClasses::Ex("There must be at least 2 folds.");
    
    model->rand().setSeed(seed);

    // Test
    std::cout.precision(8);
    double sae;
    double sse = model->repValidate(*features, *labels, reps, folds, &sae, succinct ? NULL : CrossValidateCallback, model);
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

    std::cout << " ------------ Finished Cross Validation with " << GetParameter("algorithm") << " ------------------- \n\n" ;

    return error_flag;
}

void WafflesInterface::CrossValidateCallback(void* pSupLearner, size_t nRep, size_t nFold, double foldSSE, size_t rows)
{
    std::cout << "Rep: " << nRep << ", Fold: " << nFold <<", Mean squared error: " << GClasses::to_str(foldSSE / rows) << "\n";
}

}

#endif 

