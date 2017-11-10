#include "api.h"
#include <unistd.h>

/** \file dummy.C
 * \brief Example showing the solver selecter and database builder in action for a simple problem 
 **/

static int solvetimehack = 0;
static int bad_count = 0;
static int really_bad = 0;
static int perfect_count = 0;
static std::vector< double> satime(6);
static double rtime = 0;
static double mtime = 0;

/**
 * aaaaa
 */

class CustomFeature : public _SS_Feature
{
  public:
  double datasum;

  CustomFeature() : _SS_Feature("Custom"), datasum(0.0) {}
  _SS_ErrorFlag InitializeFeatureCollection( MPI_Comm comm, const int & nrows, const int &ncols ) override 
  {
    datasum = 0;
    return _SS_error_flag;
  }
  _SS_ErrorFlag NextElement( MPI_Comm comm , const int &row, const int &col, const double &data ) override
  {
    
     datasum += data;
     return _SS_error_flag;
  } 
  virtual _SS_ErrorFlag FinalizeFeatureCollection(MPI_Comm comm, double &final_value ) override
  {
      final_value = datasum;
      return _SS_error_flag;
  }
};

class CustomMeasurement : public _SS_Measurement
{
  public:
  CustomMeasurement( const double &bad_value ) : _SS_Measurement(bad_value) 
  {
    parameter_list.insert(std::make_pair("fake_parameter",0.0134));
  }

  _SS_ErrorFlag StartMeasurement( MPI_Comm comm, std::map<std::string, double> &mstruct )
  {
      /* Set some parameters maybe */
      mstruct["fake_parameter"] = -5;
      return _SS_error_flag; /* This measurement is a hack */
  }
  _SS_ErrorFlag StopMeasurement(MPI_Comm comm, std::map<std::string, double > &mstruct, double &final_value  )
  {
     if ( mstruct["fake_parameter"] > 0 ) 
        final_value = solvetimehack;
     else
        final_value = 99999999;

     return _SS_error_flag;
  }
}; 

class DummyUI : public _SS_UserInterface<double,double>
{
  public:

  /**
   * aaaaa
   */
  DummyUI()  : _SS_UserInterface<double,double>() 
  {
  
  }

  _SS_ErrorFlag GetDataBaseImpl( std::shared_ptr< _SS_DataBaseBase > &database )
  {
    database.reset( new Database::_SS_DataBaseSql( "dummy.database" ) );
    return _SS_error_flag;
  }

  _SS_ErrorFlag GetMLImpl( std::shared_ptr< _SS_MachineLearning > &machinemodel )
  {
    machinemodel.reset( new MachineLearning::_ML_Waffles() );  
    return _SS_error_flag;
  }


   _SS_ErrorFlag GetMatrixInfo( double &A, int &nrows, int &ncols, int &nchunks, 
                                std::string &matrix_name, bool &mfree )
  {
     nrows = 1;
     ncols = 1;
     nchunks = 1;
     matrix_name = std::to_string(A);
     mfree = false;
     return _SS_error_flag;
  }

  _SS_ErrorFlag GetSparcity( double &A, int &chunk, 
                             std::vector< std::pair< unsigned int, unsigned int > > &sparcity,
                             std::vector< double > &values
                           ) override
  {
    sparcity.push_back(std::make_pair((unsigned int) 0, (unsigned int)0));
    values.push_back(A);
    return _SS_error_flag;
  }

  /** Add Features */
  _SS_ErrorFlag AddFeaturesAndMeasurements( _SS_Features &features, _SS_Measurements &measure ) override
  {
    std::shared_ptr<_SS_Feature> cf( new CustomFeature() );
    std::shared_ptr<_SS_Measurement> cm(  new CustomMeasurement(0.3) );
     features.AddFeature( cf );     
     measure.AddMeasurement("CPUTime" , cm);
     return _SS_error_flag;
  }

  
  /**
   * Solve the system as given 
   */
  _SS_ErrorFlag SolveSystem( double &A, double &x, double &b,
                             _SS_Solver &solver,
                             std::map<std::string,double> &mstruct ) override
  { 
    
    
    /* In this example, the solve time is determined by the following function 
     *
     *  ST  = ( A - c ) ^2 + c +1 if A > C 
     *
     *  where c \in { 0, 1, 2 } is the parameter "coord"  and A \in { 0,...,5 }
     *
     */

    /* Use the parameters obtained by the solver selecter */

    double coord;    
    std::string vals;
    solver.GetParameter("coord", vals);
    coord = std::stof(vals);
          
    double solvetime = ( A - coord ) * ( A - coord ) + coord + 1 ;
    double min = solvetime;
    double sstime;
    for ( int i = 0; i < 5; i++ ) 
    {       
       sstime =  ( A - i ) * ( A - i ) + i + 1 ;
       min = ( min < sstime ) ? min : sstime ;
       satime[i] += sstime;
    }
    
    /* Mstruct contains a bunch of parameters we must set. This structure
     * only needs to be set when measurements are taking place. This includes
     * inline measurements of the matrix, or when a database is being built from 
     * file. One example is that the CPUTIME should only be set if the method actually
     * converged. If it did not converge, or did not run at all, the solvetime should 
     * be set very large, but might be very short. */
    mstruct["fake_parameter"] = 100; 

    int ran = rand() % 5 ;
    sstime =  ( A - ran ) * ( A - ran ) + ran + 1 ;
    satime[5] += sstime;
    rtime += solvetime;
    mtime += min; 
    double percent = ( solvetime - min ) / min * 100; 
    if ( percent > 30 )
    { bad_count++; if ( percent  > 40 ) really_bad++; }  
    else if ( percent < 1 )
    { perfect_count++;}
    



    solvetimehack = solvetime;
    return 0;
  }

  /**
   * asdasd
   **/
  _SS_ErrorFlag InitMatrix(std::string filename, 
                           std::unique_ptr<double> &A)
  {
     /* import A from some  file */
     A.reset(new double);   *A = std::stof(filename);
     return _SS_error_flag;    
  }

  _SS_ErrorFlag InitVector( const double &A, std::unique_ptr<double> &x )
  {
      /* In this case a custom destroy function is not required. If it is, one could impliment
       * the Free*** functions in this interface.  */
      x.reset( new  double );
      *x = 0.0;
      return _SS_error_flag;
  }

  _SS_ErrorFlag SetVector( std::unique_ptr<double> &x , std::vector<int> &cols, const std::string &init_type ) override
  {
      x.reset(new double);
      if ( init_type == "random" )
          *x = (double) rand()/RAND_MAX;
      else if ( init_type == "ones" )
          *x = 1;
      else
          *x = 0;
      return _SS_error_flag;
  }
 
  
};

/** Main */
int main (int argc, char* argv[])
{

  MPI_Init(&argc,&argv);
  int build_database = 0; // should we solve a problem 
  int num_tests = 10000;
  int use_ml = 1;
  int dump = 0;
  int build_inline = 0;
  
  int i = 1;
  while ( i < argc )
  {
      if ( std::string(argv[i]) == "--buildinline" )
      {
          i++;
          build_inline = std::stoi(argv[i++]); 
      }
      else if ( std::string(argv[i]) == "--buildfromfile" )
      {
          i++;
          build_database = std::stoi(argv[i++]); 
      }
      else if ( std::string(argv[i]) == "--tests" )
      {
          i++;
          num_tests = std::stoi(argv[i++]);
      }
      else if ( std::string(argv[i]) == "--dump" )
      {
          i++;
          dump = std::stoi(argv[i++]);
      }
      else if ( std::string(argv[i]) == "--use_ml" )
      {
          i++;
          use_ml = std::stoi(argv[i++]);
      }
      else if ( std::string(argv[i]) == "--help" )
      {
          std::cout << "This is a dummy function showing the capabilities of the SolverSelecter. \n"
                    << "In this example, the \"matrix\", A, is a positive double. The solve time \n"
                    << "for A is determined by the following function \n"
                    << "\t ST  = ( A - c ) ^2 + c +1 if A > C \n"
                    << "\t where c in { 0, 1, 2,... } plays the part of a solver parameter.\n"
                    << "Two options can be specified.\n"
                    << "\t --build 1 builds the database. ( 0 solves a problem using the database ) \n"  
                    << "\t --test <int> number of solver tests to run . \n"   
                    << "\t --dump <int> dump the matrix to \"file\" . \n"   
                    << "\t --use_ml <int> turn on/off solver selection . \n"   ;
          return 0;    
      }
      else
          std::cout << "Command " << argv[i++] << " not found \n";
  }
  srand( time(NULL) );
   
  /* Make a shared pointer to the interface, and a solverselecter */
  _SS_SolverSelecter<double,double> solver_selecter;

  std::shared_ptr<DummyUI> interface( new DummyUI() );
  solver_selecter.Initialize(interface); 
  solver_selecter.SetMatrixDump(dump);
  solver_selecter.SetSolverSelection(use_ml);
  solver_selecter.SetBuildDatabaseInline(build_inline, "dummy.input"); /* test solvers on the fly */
  
  if ( build_database )
  {
    solver_selecter.BuildDataBaseFromFile("dummy.input"); /* build up data base from external list of mats */
    return 0 ;
  }

  
  else 
  {
    double A,x,b ; 
    int total = num_tests;

    if ( build_inline )
    {
        std::cout << " Building the database using inline testing based on the solver pairs defined in  \
                       the input file. If --tests <num> is to large, this might take a long time. "; 
    }
    for ( int i = 0; i < total; i++ ) 
    {
        A = (double) rand()/ (double) RAND_MAX * 5.0 ;
        solver_selecter.Solve( A, x, b );
    }
    if ( ! build_inline && !build_database )
    {
       std::cout << "Num simulations : " << total << std::endl;
       std::cout << "BAD (i.e > 30 \% of min) : " << bad_count << " ( " <<  (double) bad_count / ( total ) *100 << " \% ) \n";
       std::cout << "Really bad ( i.e > 40 \%): " << really_bad << " (" << (double) really_bad / total * 100 << " \% ) \n"; 
       std::cout << "PERFECT ( i.e == min ) : "   << perfect_count << " ( " <<  (double) perfect_count / ( total ) *100 << " \% ) \n";
       std::cout << "\n\n";    
       std::cout << "Actual Solve Time = " << rtime << std::endl;
       std::cout << "Best Possible Solve Time = " << mtime << std::endl;
       std::cout << "Percentage from best = " << (rtime - mtime)/mtime*100 << std::endl;
        
       for ( auto i = 0; i < 5; i++ )
       {
         std::cout << "Solve time using coord = " << i << " = " << satime[i] << " ( " <<  100 - (rtime/satime[i])*100 <<" \% )\n";
       }
       std::cout << "Solve time using coord = r = " << satime[5] << " ( " <<  100 - (rtime/satime[5])*100 <<" \% )\n";
    }
  }
  srand( time(NULL));

  /* The solver selecter should be finalized. This re classifies the solvers in 
   * the database when ever they are added inline, and closes the db */
  solver_selecter.Finalize();
  

  return 0;
}
  



