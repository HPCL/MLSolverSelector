
#if WITH_ZOLTAN

#include "GColoring.h"
_SS_Zoltan_Graph::_SS_Zoltan_Graph( MPI_Comm comm )
{

}


_SS_ErrorFlag _SS_Zoltan_Graph::SetOwnerShipRange( int g_num_rows, int g_num_cols)
{
    /* Global ownership ranges */
    num_global_vertices = g_num_rows + g_num_cols ;
    global_edge_list.resize(num_global_vertices);
    return _SS_error_flag;
}


_SS_ErrorFlag _SS_Zoltan_Graph::AddEdges( std::vector< std::pair< unsigned int,unsigned int > > &sparcity )
{
    for ( auto &edge : sparcity )
    {
        auto g_edge = Map2Dto1D(edge,false);
        global_edge_list[g_edge.first].push_back(g_edge.second);
        global_edge_list[g_edge.second].push_back( g_edge.first);
    }

    return _SS_error_flag;

}

void _SS_Zoltan_Graph::GetGlobalColumnID( std::vector< unsigned int > col_ids, std::vector<unsigned int> &g_ids)
{
    auto temp_pair = std::make_pair( (unsigned int) 0, (unsigned int) 0 );
    auto temp_pair1 = std::make_pair( (unsigned int) 0, (unsigned int) 0 );

    g_ids.clear();
    for (auto it : col_ids )
    {
        temp_pair.second = it;
        temp_pair1 = Map2Dto1D( temp_pair, false );
        g_ids.push_back( temp_pair1.second );
    }
}

std::pair<int,int> _SS_Zoltan_Graph::Map2Dto1D( std::pair<unsigned int,unsigned int> &inds, bool inverse)
{
    if ( !inverse )
        return std::make_pair(2*inds.first,2*inds.second+1);
    else
        return std::make_pair(inds.first/2,(inds.second-1)/2);
}

_SS_ZoltanGraphColor::_SS_ZoltanGraphColor( MPI_Comm comm , bool initialize )
{

    graph = new _SS_Zoltan_Graph( comm );

    /* This will never do anything important, but zoltan requires it to be called. All
     * MPI Init calls are handled way up the chain. Kinda annoying cause it throws a warning
     * about argv and argc being used before init. Could just init them i guess. */
    char **argv ;
    int argc;
    float ver;
    if  ( initialize )
        Zoltan_Initialize( argc, argv, &ver ); /* command line args do nothing when MPI_INIT has been called
                                             prior. */
    zz = new Zoltan( comm );
    zz->Set_Param("DEBUG_LEVEL","0");
    zz->Set_Param("NUM_GID_ENTRIES","1");
    zz->Set_Param("NUM_LID_ENTRIES","1");
    zz->Set_Param("COLORING_PROBLEM","BIPARTITE");
    zz->Set_Num_Obj_Fn( zoltan_num_obj_fn , graph );
    zz->Set_Obj_List_Fn( zoltan_obj_list_fn , graph );
    zz->Set_Num_Edges_Multi_Fn(zoltan_num_edges_multi_fn, graph );
    zz->Set_Edge_List_Multi_Fn(zoltan_edge_list_multi_fn , graph );


}

_SS_ErrorFlag _SS_ZoltanGraphColor::Clear()
{
    graph->global_edge_list.clear();
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_ZoltanGraphColor::AddEdges( std::vector< std::pair< unsigned int,unsigned int > > &sparcity )
{
    graph->AddEdges(sparcity);
    return _SS_error_flag;
}

_SS_ErrorFlag _SS_ZoltanGraphColor::ColorThatGraph( int &g_num_rows, int &g_num_cols,
        std::vector< std::pair<unsigned int,unsigned int>> &sparcity,
        std::vector< unsigned int > &column_ids,
        std::vector< int > &col_colors )
{

    graph->SetOwnerShipRange(g_num_rows,g_num_cols);

    std::vector< unsigned int > global_col_ids;
    graph->GetGlobalColumnID( column_ids , global_col_ids );

    int one(1);
    zz->Color( one, global_col_ids.size(), global_col_ids.data(), col_colors.data() );

    return _SS_error_flag;
}

int _SS_ZoltanGraphColor::zoltan_num_obj_fn( void *data, int *ierr )
{
    *ierr = ZOLTAN_OK;
    _SS_Zoltan_Graph *graph = (_SS_Zoltan_Graph*) data;
    return graph->num_global_vertices;
}

void _SS_ZoltanGraphColor::zoltan_obj_list_fn( void * data, int nge, int nle,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int wgt_dim,
        float *obj_wgts,
        int * ierr)
{
    _SS_Zoltan_Graph *graph = (_SS_Zoltan_Graph*) data;

    for ( int i = 0; i < graph->num_global_vertices; i++ )
    {
        global_ids[i] = i;
        local_ids[i] = i;
    }
    *ierr = ZOLTAN_OK;

}



 void _SS_ZoltanGraphColor::zoltan_num_edges_multi_fn( void *data, int num_gid_entries, int num_lid_entries, int num_obj,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int *num_edges,
        int *ierr )
{
    _SS_Zoltan_Graph *graph = (_SS_Zoltan_Graph*) data;

    for ( int i = 0 ; i < num_obj; i++ )
    {
        num_edges[global_ids[i]] = graph->global_edge_list[global_ids[i]].size();
    }
    *ierr = ZOLTAN_OK;
    return;

}

void _SS_ZoltanGraphColor::zoltan_edge_list_multi_fn(void *data, int num_gid_entries, int num_lid_entries, int num_obj,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int *num_edges,
        ZOLTAN_ID_PTR nbor_global_id,
        int *nbor_procs,
        int wgt_dim,
        float *ewgts,
        int *ierr )
{
    _SS_Zoltan_Graph *graph = (_SS_Zoltan_Graph*) data;
    int index = 0;

    for ( int i = 0; i < num_obj ; i++ )
    {

        for ( int j = 0; j < num_edges[i]; j++ )
        {
            nbor_global_id[index] = graph->global_edge_list[global_ids[i]][j];
            nbor_procs[index] = 0;
            index++;
        }
    }
    *ierr = ZOLTAN_OK;
    return;

}

_SS_ZoltanGraphColor::~_SS_ZoltanGraphColor()
{
    delete graph;
    delete zz;
}

#endif
