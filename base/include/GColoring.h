#ifndef GCOLORING_H
#define GCOLORING_H

#if WITH_ZOLTAN

#include "typedefs.h"
#include "zoltan.h"
#include "zoltan_cpp.h"

class _SS_Zoltan_Graph
{
public:
    int num_global_vertices;
    std::vector< std::vector< unsigned int > > global_edge_list;

    _SS_Zoltan_Graph( MPI_Comm comm );

    _SS_ErrorFlag SetOwnerShipRange( int g_num_rows, int g_num_cols);

    _SS_ErrorFlag AddEdges( std::vector< std::pair< unsigned int,unsigned int > > &sparcity );

    void GetGlobalColumnID( std::vector< unsigned int > col_ids, std::vector<unsigned int> &g_ids);

    std::pair<int,int> Map2Dto1D( std::pair<unsigned int,unsigned int> &inds, bool inverse);

};

class _SS_ZoltanGraphColor
{
public:
    _SS_Zoltan_Graph *graph;
    Zoltan* zz;
    int rank,size;

    _SS_ZoltanGraphColor( MPI_Comm comm , bool initialize );

    ~_SS_ZoltanGraphColor();

    _SS_ErrorFlag Clear() ;

    _SS_ErrorFlag AddEdges( std::vector< std::pair< unsigned int,unsigned int > > &sparcity );

    _SS_ErrorFlag ColorThatGraph( int &g_num_rows, int &g_num_cols,
                                  std::vector< std::pair<unsigned int,unsigned int>> &sparcity,
                                  std::vector< unsigned int > &column_ids,
                                  std::vector< int > &col_colors );

    static int zoltan_num_obj_fn( void *data, int *ierr );

    static void zoltan_obj_list_fn( void * data, int nge, int nle,
                                    ZOLTAN_ID_PTR global_ids,
                                    ZOLTAN_ID_PTR local_ids,
                                    int wgt_dim,
                                    float *obj_wgts,
                                    int * ierr) ;


    static void zoltan_num_edges_multi_fn( void *data, int num_gid_entries, int num_lid_entries, int num_obj,
                                           ZOLTAN_ID_PTR global_ids,
                                           ZOLTAN_ID_PTR local_ids,
                                           int *num_edges,
                                           int *ierr );

    static void zoltan_edge_list_multi_fn(void *data, int num_gid_entries, int num_lid_entries, int num_obj,
                                          ZOLTAN_ID_PTR global_ids,
                                          ZOLTAN_ID_PTR local_ids,
                                          int *num_edges,
                                          ZOLTAN_ID_PTR nbor_global_id,
                                          int *nbor_procs,
                                          int wgt_dim,
                                          float *ewgts,
                                          int *ierr );

};

#endif
#endif
