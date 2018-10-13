#ifndef MOOSEINTERFACE_H
#define MOOSEINTERFACE_H

#ifdef WITH_PETSCUI

#include "PetscInterface.h"

namespace SolverSelecter {


class MooseCoupler : public PetscCoupler {
  public:
  static bool commandLine;
  static std::string inputFile, serialized_model, algorithm, database;
  static std::string useML, dumpMatrix, buildInline; 

  static std::map< std::string, std::map< std::string, std::set< std::string > > > solvers;
  static std::map< std::string, std::map< std::string, std::set< std::string > > > precons;
  static std::set< std::pair< std::string, std::string > > solver_pairs;
  
  static void setupCoupler(  std::string _input, std::string _data, std::string _serial, std::string _alg,  std::string _use, std::string _dump, std::string _build );


  static void addSolver( std::string name );
  static void addPreconditioner( std::string name ); 
  static void addSolverParameter( std::string solver, std::string parameter, std::vector< std::string > &values ); 
  static void addPreconditionerParameter( std::string solver, std::string parameter, std::vector< std::string > &values ); 
  static void addParameter( std::map<std::string, std::map<std::string, std::set<std::string>> > &m, 
                            std::string name, std::string parameter, std::vector< std::string > &values );
  static void addSolverPair( std::string solver, std::string precon );  
  
  static PetscErrorCode KSPCreate_MooseSS(KSP ksp);
  static PetscErrorCode KSPSetFromOptions_MooseSS(PetscOptionItems *PetscOptionsObject,KSP ksp); 
  static void RegisterKSP();


    
private:    
    
    static void addM( std::map<std::string, std::map< std::string, std::set< std::string >>> &m, std::string name );     
    static void addP( std::map<std::string, std::map<std::string, std::set< std::string >>> &m ,  
                      std::string s, 
                      std::string s1, 
                      std::vector<std::string> s2 ) ;

};

class MooseFileParser : public AsciiFileParser {

  ErrorFlag Parse(std::vector< std::string > &filenames , std::vector< Solver > &solver_list ) override;
};

class MooseInterface : public PetscUI
{
public:
 
    virtual ErrorFlag GetInputFileImpl( std::shared_ptr< InputFileInterface > &parser ) override;    
  
};

}


#endif 
#endif

