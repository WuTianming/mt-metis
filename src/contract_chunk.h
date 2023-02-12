/**
 * @file contract_swap.c
 * @brief Functions for performing contraction with disk-based adjlist chunks
 * @author wtm
 * @version 1
 */




#ifndef MTMETIS_CONTRACT_CHUNK_H
#define MTMETIS_CONTRACT_CHUNK_H



#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


void par_contract_chunk_graph(
    ctrl_type * ctrl, 
    graph_type * graph, 
    vtx_type mycnvtxs, 
    vtx_type const * const * gmatch, 
    vtx_type const * fcmap);




#endif
