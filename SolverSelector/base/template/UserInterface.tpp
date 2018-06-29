
namespace SolverSelecter
{

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::ClassificationMeasurements( std::map<std::string, double> &class_values )
{
    std::cout << "Classification Measurements not implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::DumpMatrix(Matrix &A )
{
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::InitMatrix(std::string filename,
        std::unique_ptr<Matrix> &A )
{
    std::cout << "Init Matrix Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::InitVector(const Matrix &AA,
        std::unique_ptr<Vector> &x )
{
    std::cout << "Init Vector Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::CopyVector(const Vector &AA,
        std::unique_ptr<Vector> &x )
{
    std::cout << "Copy Vector Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::SetVector(std::unique_ptr<Vector> &cloned_x, std::string type_ )
{
    /* type is either "ones" for a vector of ones, or "rand" for a random vector. */
    std::cout << "SetVector Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag UserInterface<Matrix,Vector>::FreeMatrix( std::unique_ptr<Matrix> &A )
{
    std::cout << "Free Matrix Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::FreeVector( std::unique_ptr<Vector> &x )
{
    std::cout << "Free Vector Not Implimented\n";
    return error_flag;
}

template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::GetMatrixName( Matrix &A , std::string &name )
{
    name = "default_name" + std::to_string( std::rand() );
    return error_flag;
}


template<typename Matrix, typename Vector>
ErrorFlag
UserInterface<Matrix,Vector>::GetDefaultSolver( Solver &solver )
{
    std::cout << "GetDefaultSolver Not Implimented\n";
    return error_flag;
}

}
