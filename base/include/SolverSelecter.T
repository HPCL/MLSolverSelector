

/** \file SolverSelecter.t
 * \brief Templete implimentation file for the solver selecter class
 */
template <typename Matrix, typename Vector>
_SS_SolverSelecter<Matrix,Vector>::_SS_SolverSelecter() : _measurements(MPI_COMM_WORLD), _features(MPI_COMM_WORLD)
{
    Add("use_ml","true","Should we use ML to determine best SS. Otherwise, act as a simple pass through");
    Add("matrix_dump","false","Dump the matrix to file prior to solving");
    Add("DBInline","false","Do inline measurements for database building");
    Add("matrix_free","false","use matrix free extraction techniques");

#if WITH_ZOLTAN
    graph = new _SS_ZoltanGraphColor(MPI_COMM_WORLD,true);
#endif
}


template <typename Matrix, typename Vector>
_SS_SolverSelecter<Matrix,Vector>::~_SS_SolverSelecter()
{
#if WITH_ZOLTAN
    delete graph;
#endif
}

/**  Initialize the solver selecter.  */
template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::Initialize(std::shared_ptr<_SS_UserInterface<Matrix,Vector>> _interface)
{
    interface = _interface;
    interface->SetSolverSelecterOptions(this);
    interface->AddFeaturesAndMeasurements(_features,_measurements);
    interface->GetDataBaseImpl( database );
    interface->GetMLImpl( machinemodel );

    _measurements.GetParameterList(_mstruct);
    database->Initialize();
    machinemodel->Initialize(database);

    return _SS_error_flag;
}


/** Finalize the solver selecter  */
template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::Finalize()
{
    std::string DBInline;
    Get("DBInline",DBInline);
    if ( DBInline == "true" )
        ClassifySolvers();

    database->Finalize();
    return _SS_error_flag;
}


template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::SetSolverSelection(bool on)
{
    if (on)
        Set("use_ml","true");
    else
        Set("use_ml","false");

    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::SetMatrixDump(bool on)
{
    if (on)
        Set("matrix_dump","true");
    else
        Set("matrix_dump","false");

    return _SS_error_flag;
}



/** Get features for a matrix free matrix */
template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::ExtractFeatures( Matrix &A, std::string &mfile )
{

    sparcity.clear();
    values.clear();

    int sparcity_chunks;
    bool matrix_free;
    interface->GetMatrixInfo(A, nrows, ncols, sparcity_chunks, mfile, matrix_free );

    ones.resize(nrows);/*TODO What does this do, AND WHY IS IT HERE */

    if ( !matrix_free )
    {
        _features.Initialize(nrows,ncols);
        for ( int i = 0; i < sparcity_chunks ; i++ )
        {
            interface->GetSparcity( A, i, sparcity, values );
            for ( unsigned int i = 0; i < sparcity.size(); i++ )
            {
                _features.Next( sparcity[i].first, sparcity[i].second, values[i] );
            }
            sparcity.clear();
            values.clear();
        }
        _features.Finalize();
    }
    else
    {
#if WITH_ZOLTAN

        graph->Clear();
        graph->graph->SetOwnerShipRange(nrows,ncols);
        for ( int i = 0; i < sparcity_chunks; i++ )
        {
            interface->GetSparcity( A, i, sparcity, values );
            graph->AddEdges(sparcity);
            sparcity.clear();
            values.clear();
        }

        std::vector<unsigned int> columns_needed(ncols);
        std::vector<int> column_cols(ncols);
        std::iota( std::begin(columns_needed), std::end(columns_needed), 0);
        graph->ColorThatGraph( nrows, ncols, sparcity, columns_needed, column_cols );
        graph->Clear();

        std::vector< std::vector< int > > col_vecs;
        for ( unsigned int i = 0 ; i < column_cols.size() ; i++ )
        {
            if (col_vecs.size() < (unsigned int) column_cols[i] )
                col_vecs.resize( column_cols[i] );
            if ( col_vecs[column_cols[i]-1].size() != (unsigned int) nrows )
                col_vecs[column_cols[i]-1].resize(nrows,0);
            col_vecs[column_cols[i]-1][i] = 1;
        }

        interface->InitVector( A, b );

        std::vector<double> vals(nrows);
        std::vector<int> rows(nrows);
        std::iota( std::begin(rows), std::end(rows) , 0 );

        _features.Initialize(nrows,ncols);

        for ( auto it : col_vecs )
        {

            interface->SetVector( b, it, "ones" );
            interface->MatVecAndGetData(A,*b, rows, vals );
            for ( unsigned int i = 0; i < rows.size(); i++ )
            {
                _features.Next( rows[i], -1 , vals[i] );
            }
        }
        _features.Finalize();

        interface->FreeVector( b );
#else
        std::cout << " Matrix free not supported without Zoltan support. Who knows what will happen now??? \n ";
#endif

    }
    return _SS_error_flag;

}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::SetBuildDatabaseInline( bool on,  std::string input_file )
{
    if ( on )
    {
        ParameterSpaceFromFile(input_file); /* Gets the solvers from the input file */
        Set("DBInline", "true"); /* Set the flag to do inline DB measurement */
    }
    else
        Set("DBInline", "false");

    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::BuildDataBaseFromFile( const std::string &inputfile )
{
    ParameterSpaceFromFile( inputfile );
    for ( auto matrix_file : filenames )
    {
        interface->InitMatrix( matrix_file, A);
        interface->InitVector( *A, x );
        interface->InitVector( *A, b );

        ExtractFeatures( *A, matrix_file );
        _mfile = matrix_file;

        for ( unsigned int i = 0; i < solvers.size(); i++ )
        {
            _solver.Clear();
            _solver.SetSolverName(solvers[i],preconds[i]);
            for (auto param : parameters[i] )
            {

                std::string value;
                if ( param.second.size() > 0 )
                    value =  *( (param.second).begin() );
                _solver.SetParameter( param.first, value );
            }

            //Loop over parameter space
            RecurseParameterSpace( *A, parameters[i] );
        }

        interface->FreeMatrix(A);
        A.reset(nullptr);
        interface->FreeVector(x);
        x.reset(nullptr);
        interface->FreeVector(b);
        b.reset(nullptr);
    }

    ClassifySolvers();
    database->WriteToFile();
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::GetDataBase( _SS_DataBaseBase **_dat )
{
    *_dat = database;
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
std::shared_ptr<_SS_UserInterface<Matrix,Vector>> 
_SS_SolverSelecter<Matrix,Vector>::GetInterface()
{
    return interface;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::GetMachineLearner( _SS_MachineLearning *machine ) const
{
    machine = &machinemodel;
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::AddAllowedSolver( const std::string &solver, const std::string &precond )
{
    sset.push_back( std::make_pair( solver, precond ) );
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::PrintChosenSolver()
{
    _solver.Print();
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::Solve( Matrix &AA, Vector &xx, Vector &bb )
{
    std::string use_ml, dump_mat, DBInline;
    Get("use_ml",use_ml);
    Get("DBInline",DBInline);

    Get("matrix_dump",dump_mat);
    if ( dump_mat == "true" )
    {
        interface->DumpMatrix( AA );
    }

    if ( use_ml == "true" && DBInline == "false" )
    {
        /* Train the model based on the data base (if not done already) */ 
        auto start = std::chrono::high_resolution_clock::now();
        machinemodel->Train( sset );
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        std::cout << " Trained the system in " << (double) microseconds/1000000 << " seconds \n" ;

        start = std::chrono::high_resolution_clock::now();
        ExtractFeatures( AA, _mfile );
        elapsed = std::chrono::high_resolution_clock::now() - start;
        microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        std::cout << " Extracted Features in " << (double) microseconds/1000000 << " seconds \n" ;
        

        /* Classify the matrix and get the solver */
        start = std::chrono::high_resolution_clock::now();
        machinemodel->Classify(_features, _solver );
        elapsed = std::chrono::high_resolution_clock::now() - start;
        microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        std::cout << " Classified the system in " << (double) microseconds/1000000 << " seconds \n" ;
        

        start = std::chrono::high_resolution_clock::now();
        interface->SolveSystem( AA, xx, bb, _solver, _mstruct );
        elapsed = std::chrono::high_resolution_clock::now() - start;
        microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        std::cout << " Solved the system in " << (double) microseconds/1000000 << " seconds \n" ;
        

    }
    else if ( DBInline == "true" )
    {
        /* Extract the Features from the matrix. ( user implimented ) */
        ExtractFeatures( AA,  _mfile );
        interface->InitVector( AA, x );
        interface->InitVector( AA, b );

        for ( unsigned int i = 0; i < solvers.size(); i++ )
        {
            _solver.Clear();
            _solver.SetSolverName(solvers[i],preconds[i]);
            for (auto param : parameters[i] )
            {

                std::string value;
                if ( param.second.size() > 0 )
                    value =  *( (param.second).begin() );
                _solver.SetParameter( param.first, value );
            }
            //Loop over parameter space
            RecurseParameterSpace( AA, parameters[i] );
        }
        interface->FreeVector(x);
        x.reset(nullptr);
        interface->FreeVector(b);
        b.reset(nullptr);

        /*TODO Currently, just use the last solver to solve the actual problem. Best approach
         * would be to use machine learning to solve the actaual problem. This means we need to
         * somehow classify the solvers above as good or bad for this matrix problem.
         */
        interface->SolveSystem(AA,xx,bb,_solver, _mstruct );

    }
    else
    {   
        std::vector< int > hash_list;
        database->GetSolverHashList(hash_list);
        database->GetSolver( hash_list[0] , _solver );       
        interface->SolveSystem(AA,xx,bb,_solver, _mstruct);
    }

    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::Train()
{

    std::string use_ml, dump_mat;
    Get("use_ml",use_ml);
    /* Only train if ML is turned on, else, do nothing */
    if ( use_ml == "true" )
    {
        machinemodel->Train( sset );
    }
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::RecurseParameterSpace(Matrix &AA,  std::map< std::string, std::set< std::string >> &solver_params )
{
    if ( solver_params.size() == 0 )
    {
        interface->SetVector( x, ones, "random" );
        interface->SetVector( b, ones, "ones" );

        _measurements.Start(_mstruct);
        interface->SolveSystem(AA,*x,*b,_solver, _mstruct );
        _measurements.Stop(_mstruct);

        WriteToDataBase();

        return _SS_error_flag;
    }

    auto its = std::prev( solver_params.end() );
    if ( its->second.size() > 0 )
    {
        for ( auto it : its->second )
        {
            _solver.SetParameter( its->first, it );
            auto sparams = solver_params;
            sparams.erase(std::prev( sparams.end() ) );
            RecurseParameterSpace(AA, sparams);
        }
    }
    else
    {
        auto sparams = solver_params;
        sparams.erase(std::prev(sparams.end()));
        RecurseParameterSpace(AA, sparams);
    }
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::WriteToDataBase( )
{
    database->AddRow( _solver, _mfile, _features, _measurements );
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::ClassifySolvers( )
{
    /* these will be parameters, 1 indicates twice as fast as
    * the fastest is ok with us */
    database->ClassifySolvers( _measurements );

    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::ParameterSpaceFromFile(std::string fname )
{

    std::fstream file;
    file.open(fname.c_str(), std::ios::in | std::ios::out | std::ios::app | std::ios::binary );

    std::string line, l, hypre_subname;
    std::vector <std::string > split_line, result;

    //adding = -1 for something else, 0 for solver, 1 for precond, 2 for hypre,
    int adding = -1;

    std::string solver_name;
    std::map< std::string, std::set< std::string > > solver_params;

    std::map< std::string, std::map< std::string, std::set< std::string > > > solvers;
    std::map< std::string, std::map< std::string, std::set< std::string > > > precons;
    std::map< std::string, std::map< std::string, std::set< std::string > > > hypres;
    std::set< std::pair< std::string, std::string > > solver_pairs;

    while (!file.eof())
    {
        std::getline(file,line); /* get the next line in the file */
        result.clear();

        StringSplit( line, " ", result );

        if ( result.size() == 0 )
        {
            /* Do nothing */
        }
        else if ( result[0] == "@PARAMETER" )
        {
            if ( adding < 0 )
                std::cout << " File in wrong format. Adding parameter without solver \n" ;
            else
            {
                std::string pname = result[1];
                std::set< std::string > pset;
                if ( result.size() <= 2 )
                    pset.insert("none");
                else
                {
                    for ( unsigned int i = 2; i < result.size(); i++ )
                    {
                        std::string fi(result[i].begin(),result[i].begin()+1);
                        if ( fi == "#" ) break;
                        pset.insert(result[i]);
                    }
                }
                solver_params.insert( std::make_pair(pname, pset) );
            }
        }
        else
        {
            /* Finish up the previous solver if required */
            if ( adding == 0 )
                solvers.insert(std::make_pair( solver_name, solver_params ) );
            else if ( adding == 1 )
                precons.insert(std::make_pair( solver_name, solver_params ) );
            else if ( adding == 2 )
                hypres.insert(std::make_pair( hypre_subname, solver_params ) );
            adding = -1;

            if ( result[0] == "@MATRIX" )
            {
                AddMatrix(result[1]);
            }
            else if ( result[0] == "@PAIR" )
            {
                for ( unsigned int i = 2; i < result.size(); i++ )
                    solver_pairs.insert( std::make_pair( result[1],result[i] ) );
            }
            else if ( result[0] == "@SOLVER" )
            {
                solver_name = result[1];
                solver_params.clear();
                adding = 0;
            }
            else if ( result[0] == "@PRECON" )
            {
                solver_name = result[1];
                solver_params.clear();
                if ( solver_name == "hypre" )
                {
                    adding = 2;
                    hypre_subname = result[2];
                }
                else
                    adding = 1;
            }
            else
            {
                std::string fi(result[0].begin(),result[0].begin()+1);
                if ( fi != "#" )
                    std::cout << " Keyword " << result[0] << " not recongnized \n";
            }
        }
    }
    /* Finally, add all the solvers */
    for ( auto it : solver_pairs )
    {
        if ( it.second == "hypre" )
            std::cout << " Hypre solvers refurred to by subname in solver pairs \n ";

        auto search_solver = solvers.find(it.first);
        auto search_precon = precons.find(it.second);

        if ( search_solver != solvers.end() )
        {

            if ( search_precon != precons.end() )
                AddSolver( it.first, it.second, search_solver->second, search_precon->second );
            else
            {
                auto search_hypre = hypres.find(it.second);

                if ( search_hypre != hypres.end() )
                {
                    (search_hypre->second)["-pc_hypre_type"].insert(search_hypre->first);
                    AddSolver( it.first, "hypre", search_solver->second, search_hypre->second );
                }
                else
                    std::cout << " precon " << it.second << " must be added to input file \n ";
            }
        }
        else
            std::cout << " solver " << it.first << " must be added to input file \n ";
    }

    /* Last, do a check to make sure we have at least one solver and one matrix. */
    if ( filenames.size() == 0 || solvers.size() == 0  ) 
    {
        std::cout << " Does the file contain a solver and a matrix ???. Also, Error handling? TODO " << std::endl;
    }
    return _SS_error_flag;
}


template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::AddMatrix( std::string filename )
{
    filenames.insert(filename);
    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::AddSolver( std::string solver,
        std::string preconditioner,
        std::map< std::string, std::set< std::string >> &sparameters,
        std::map< std::string, std::set< std::string >> &pparameters )
{
    solvers.push_back( solver );
    preconds.push_back( preconditioner );

    std::map< std::string, std::set< std::string >> local;
    parameters.push_back( local );
    for ( auto it : sparameters )
        parameters.back().insert(it);
    for ( auto it : pparameters )
        parameters.back().insert(it);

    return _SS_error_flag;
}

template <typename Matrix, typename Vector>
_SS_ErrorFlag _SS_SolverSelecter<Matrix,Vector>::StringSplit( const std::string &s, const char *delim, std::vector< std::string > &result )
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while ( std::getline( ss,item,delim[0]) )
    {
        if (!item.empty())
            result.push_back( item );
    }
    return _SS_error_flag;
}

