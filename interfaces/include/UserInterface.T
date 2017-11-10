
/** \file UserInterface.T
 * \brief Templete implimentation file for the user interface class
 */

template<typename Matrix,typename Vector>
_SS_UserInterface<Matrix,Vector>::_SS_UserInterface() {};

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::AddFeaturesAndMeasurements( _SS_Features &features, _SS_Measurements &measure )
{
    std::cout << " Please override AddFeatures and add some features and measurements\n" ;
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::SetSolverSelecterOptions( _SS_SolverSelecter<Matrix,Vector> *ss )
{
    return _SS_error_flag;
}


template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::GetMatrixInfo( Matrix &A, int &nrows, int &ncols, int &chunks, std::string &matrix_name, bool &mfree )
{
    std::cout << " Impliment GetMatrixInfo please. nrows is the global number of rows, ncols is the "
              << " global number of columns, sparcity is the indecies (global) with non-zeros, and "
              << " values is the values of those non-zeros (must be in same order). If values is left "
              << " empty, then the zoltan graph coloring algorithm is used to fill in the values (for "
              << " matrix free methods). In that case, MatVec must also be implimented.  \n " ;
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::GetSparcity( Matrix &A, int &chunk,
        std::vector< std::pair< unsigned int, unsigned int > > &sparcity,
        std::vector< double > &values
                                            )
{
    std::cout << " Impliment GetSparcity \n ";
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::MatVecAndGetData( Matrix &A, Vector &x, std::vector<int> &rows, std::vector<double> &vals )
{
    std::cout << "Impliment MatVec or set values in GetMatrixInfo \n";
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::DumpMatrix( Matrix &A )
{
    std::cout << " Impliment DumpMatrix to use this feature \n";
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::_SS_UserInterface::InitMatrix( std::string filename,
        std::unique_ptr<Matrix> &A
                                                              )
{
    std::cout << " Impliment InitSystem to use the Build from file feature. \n";
    return _SS_error_flag;
}


template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::InitVector( const Matrix &AA,
        std::unique_ptr<Vector> &x )
{
    std::cout << " Impliment init vector to use the database building features \n";
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::FreeMatrix( std::unique_ptr<Matrix> &A )
{
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::FreeVector( std::unique_ptr<Vector> &x )
{
    return _SS_error_flag;
}

template<typename Matrix,typename Vector>
_SS_ErrorFlag _SS_UserInterface<Matrix,Vector>::SetVector( std::unique_ptr<Vector> &cloned_x, std::vector<int> &cols, const std::string &type )
{
    std::cout << " Impliment SetVector to use the Build database features \n";
    return _SS_error_flag;
}

