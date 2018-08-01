

#ifdef WITH_PETSCUI
#include "MooseInterface.h"

namespace SolverSelecter
{
  bool MooseCoupler::commandLine = false; 
  std::string MooseCoupler::database;
  std::string MooseCoupler::inputFile;
  std::string MooseCoupler::useML;
  std::string MooseCoupler::dumpMatrix;
  std::string MooseCoupler::buildInline;
  std::string MooseCoupler::serialized_model;
  std::string MooseCoupler::algorithm;  
  std::map< std::string, std::map< std::string, std::set< std::string > > > MooseCoupler::solvers;
  std::map< std::string, std::map< std::string, std::set< std::string > > > MooseCoupler::precons;
  std::set< std::pair< std::string, std::string > > MooseCoupler::solver_pairs;


  void MooseCoupler::setupCoupler( std::string _input, std::string _data, std::string _serial, std::string _alg,  std::string _use, std::string _dump, std::string _build ) {
        inputFile = _input;
        serialized_model = _serial;
        algorithm = _alg;
        useML = _use;
        dumpMatrix = _dump;
        buildInline = _build;
        database = _data;
        
    }

  void MooseCoupler::addSolverParameter( std::string name, std::string parameter, std::vector<std::string> &values ) 
  {
      if ( solvers.find( name ) != solvers.end() ) {
         addParameter( solvers , name, parameter, values ) ; 
      }  else {
         std::cout << "Error Solver or preconditioner does not exist" << std::endl;
      } 
  }   
  
  void MooseCoupler::addPreconditionerParameter( std::string name, std::string parameter, std::vector<std::string> &values ) 
  {
      if ( precons.find( name ) != precons.end() ) {
         addParameter( precons , name, parameter, values ) ; 
      }  else {
         std::cout << "Error preconditioner does not exist" << std::endl;
      } 
  }   
  
  void MooseCoupler::addParameter( std::map< std::string, std::map< std::string, std::set< std::string > > > &m,
                                   std::string name, 
                                   std::string parameter,
                                   std::vector< std::string> &values) {
      
      if ( m.find(name) == m.end() ) return ; 
      
      if ( m[name].find(parameter) == m[name].end() ) {
          m[name].insert( std::make_pair( parameter, std::set<std::string>() ));
      }
      for ( auto &it : values ) {
          m[name][parameter].insert(it);
      }
  } 

  void MooseCoupler::addSolver( std::string name ) {
      if ( solvers.find( name ) == solvers.end() ) {  
        solvers.insert( std::make_pair( name, std::map<std::string, std::set<std::string>>() ));
      }
  }
  
  void MooseCoupler::addPreconditioner( std::string name ) {
    if ( precons.find( name ) == precons.end() ) {
      precons.insert( std::make_pair( name, std::map<std::string, std::set<std::string>>() ));
    }
  }

  void MooseCoupler::addSolverPair( std::string solver, std::string precon ) {
        solver_pairs.insert( std::make_pair( solver, precon ) ); 
  } 
     
   PetscErrorCode MooseCoupler::KSPCreate_MooseSS(KSP ksp)
  {
      PetscFunctionBegin;
      ksp->ops->setup          = KSPSetUp_SS;
      ksp->ops->solve          = KSPSolve_SS;
      ksp->ops->destroy        = KSPDestroy_SS;
      ksp->ops->view           = KSPView_SS;
      ksp->ops->setfromoptions = KSPSetFromOptions_MooseSS;
      ksp->ops->buildsolution  = KSPBuildSolutionDefault;
      ksp->ops->buildresidual  = KSPBuildResidualDefault;

      KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,4);
      KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_RIGHT,3);
      KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_SYMMETRIC,2);
      KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_RIGHT,1);
      KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_LEFT,1);

      PetscSS *ss = new PetscSS( std::make_shared< MooseInterface >() );
      ksp->data = (void*) ss;

      PetscFunctionReturn(0);
  }


  PetscErrorCode MooseCoupler::KSPSetFromOptions_MooseSS(PetscOptionItems *PetscOptionsObject,KSP ksp)
  {
    PetscFunctionBegin;
    PetscSS *ss = (PetscSS*) ksp->data ;
    PetscUI *ui = (PetscUI*) ss->interface.get();
    bool flag = false;

    PetscOptionsHead(PetscOptionsObject,"KSP SS options");

    PetscBool flg;

    char inname[256];
    PetscStrcpy(inname, ui->input_file.c_str());
    PetscOptionsString("-ksp_ss_inputfile", " Name of the inputfile", NULL, inname,inname,256,&flg);
    if (flg) {
      ui->input_file = inname;
      MooseCoupler::commandLine = true;
      ss->Initialize(ui->input_file);
    } else {
      ui-> input_file = "Not_set_on_the_command_line_i_hope_its_in_the_input";
      ss->Initialize(ui->input_file);
    } 

    PetscOptionsTail();
    PetscFunctionReturn(0);
  }

   void MooseCoupler::RegisterKSP() {
          KSPRegister("KSPSS", KSPCreate_MooseSS);
  }

  ErrorFlag
  MooseFileParser::Parse(std::vector< std::string > &filenames ,
                       std::vector< Solver > &solver_list )
  {

      //File names is ignored for the moose parser. Bascially, 

      if ( MooseCoupler::commandLine ) { 
          AsciiFileParser::Parse(filenames, solver_list);
      } else if ( ! MooseCoupler::inputFile.empty() ) {
        inputfile = MooseCoupler::inputFile ;
        AsciiFileParser::Parse(filenames, solver_list);
      } else {
        
        std::cout << "sdfsdfsdf " << MooseCoupler::solver_pairs.size() << std::endl;
        
        GetAllSolvers(MooseCoupler::solver_pairs, MooseCoupler::solvers, MooseCoupler::precons,solver_list );

        database->SetParameter( "Database", MooseCoupler::database );
        interface->SetParameter( "build_inline", MooseCoupler::buildInline );
        interface->SetParameter( "dump_matrix", MooseCoupler::dumpMatrix );
        interface->SetParameter( "use_ml", MooseCoupler::useML );
        machinemodel->SetParameter("algorithm", MooseCoupler::algorithm);
        machinemodel->SetParameter("serialized_model", MooseCoupler::serialized_model);        
      }
      return error_flag;
  }


  ErrorFlag MooseInterface::GetInputFileImpl( std::shared_ptr< InputFileInterface > &parser ) {
      parser.reset( new MooseFileParser() );
      return error_flag;

 } 


}

#endif
