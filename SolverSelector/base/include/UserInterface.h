#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Solvers.h"
#include "MachineLearningInterface.h"
#include "DatabaseInterface.h"
#include "InputFileInterface.h"

namespace SolverSelecter { 

template<typename Matrix,typename Vector>
class UserInterface : public SSBase
{
public:
   
    
    UserInterface() : SSBase("UserInterface") 
    {
    
      AddParameter("matrix_dump","false","Dump matrix");
      AddParameter("use_ml","true","Use MachineLearning");
      AddParameter("build_inline","false","Build Database inline");
    
    } ;

    bool GetDumpMatrix() { return GetParameterAsBool("matrix_dump") ; } 
    bool GetUseML() { return GetParameterAsBool("use_ml"); }
    bool GetBuildInline() { return GetParameterAsBool("build_inline"); }
    

    virtual ErrorFlag
    SolveSystem(Matrix &A,
                Vector &x, 
                Vector &b, 
                Solver &solver, 
                bool &success ) = 0 ;

    
    virtual ErrorFlag
    ChangeSolver(Matrix &A,
                 Vector &x,
                 Vector &b ,
                 bool &change_it ) = 0;

    virtual ErrorFlag
    ExtractFeatures(Matrix &A,
                    std::map<std::string, double> &fmap ) = 0;
   
    virtual ErrorFlag
    StartMeasurements(Matrix &A,
                      Vector &x,
                      Vector &b ) = 0;

    virtual ErrorFlag
    StopMeasurements(Matrix &A,
                     Vector &x,
                     Vector &b,
                     std::map<std::string, double> &mmap ) = 0;
    
    virtual ErrorFlag
    GetDataBaseImplimentation(std::shared_ptr< DatabaseInterface > &database) = 0;
    
    virtual ErrorFlag
    GetMachineLearningModel(std::shared_ptr< MachineLearningInterface > &machine) = 0;

    virtual ErrorFlag
    GetInputFileImpl( std::shared_ptr< InputFileInterface > &parser ) = 0;


    virtual ErrorFlag
    ClassificationMeasurements( std::map<std::string, double> &class_values ) ; 
    
    virtual ErrorFlag
    DumpMatrix(Matrix &A ) ;

    virtual ErrorFlag
    InitMatrix(std::string filename,
               std::unique_ptr<Matrix> &A );

    virtual ErrorFlag
    InitVector(const Matrix &AA,
               std::unique_ptr<Vector> &x );
    
    virtual ErrorFlag
    CopyVector(const Vector &AA,
               std::unique_ptr<Vector> &x );
    
    virtual ErrorFlag
    SetVector(std::unique_ptr<Vector> &cloned_x, std::string type_ ) ;

    virtual ErrorFlag FreeMatrix( std::unique_ptr<Matrix> &A ) ;

    virtual ErrorFlag
    FreeVector( std::unique_ptr<Vector> &x ) ;
    
    virtual ErrorFlag
    GetMatrixName( Matrix &A , std::string &name );

    virtual ErrorFlag
    GetDefaultSolver( Solver &solver ); 

};

}

#include "UserInterface.tpp"

#endif
