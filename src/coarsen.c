/**
 * @file coarsen.c
 * @brief Coarsening functions
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-03
 */




#ifndef MTMETIS_COARSEN_C
#define MTMETIS_COARSEN_C




#include "coarsen.h"
#include "aggregate.h"
#include "contract.h"
#include "check.h"




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


graph_type * par_coarsen_graph(
    ctrl_type * ctrl,
    graph_type * graph)
{
  vtx_type cnvtxs;
  vtx_type ** gmatch;
  vtx_type * fcmap;
  vtx_type * match;

  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);
  tid_type const myid = dlthread_get_id(ctrl->comm);

  if (graph->nvtxs <= 1) {
    return graph;
  }

  DL_ASSERT_EQUALS((size_t)ctrl->comm,(size_t)graph->comm,"%zu");

  if (myid == 0) { 
    dl_start_timer(&ctrl->timers.coarsening);
  }

  fcmap = vtx_alloc(graph->mynvtxs[myid]);
  match = vtx_init_alloc(NULL_VTX,graph->mynvtxs[myid]);

  // DL_ASSERT(check_graph(graph),"Invalid graph");

  gmatch = dlthread_get_shmem(sizeof(vtx_type*)*nthreads,ctrl->comm);

  gmatch[myid] = match;   // distributed initialization

  par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"Graph{%zu} has %" \
      PF_VTX_T"[%"PF_VTX_T"] vertices, %"PF_ADJ_T" edges, and %"PF_TWGT_T \
      " exposed edge weight.\n",graph->level,graph->nvtxs,ctrl->coarsen_to, \
      graph->nedges,graph->tadjwgt);

  /* allocate memory for cmap, if it has not already been done due to
     multiple cuts */
  if (myid == 0) {
    if (graph->cmap == NULL) {
      /* threads need to allocate their own chunk inside the matching 
       * functions */
      graph->cmap = r_vtx_alloc(nthreads);
    }

    /* set the maximum allowed coarsest vertex weight */
    for (int i = 0; i < graph->ncon; ++i) {
      ctrl->maxvwgt[i] = \
          /* 1.5 * graph->tvwgt / ctrl->coarsen_to; */
          // 1.5*graph->tvwgt[i] / dl_max(ctrl->coarsen_to,graph->nvtxs/4.0);

          // when dealing with training_mask, tvwgt/nvtxs * 4.0 is too small (<1)
          // so the aggregation step will just be skipped
          // 1.5*(graph->tvwgt[i] / dl_max(ctrl->coarsen_to,graph->nvtxs/4.0) + 2);
          1.5*(graph->tvwgt[i] / (ctrl->coarsen_to) + 2);
    }
  }
  dlthread_barrier(ctrl->comm);

  cnvtxs = par_aggregate_graph(ctrl,graph,gmatch,fcmap);    // find a match

  int single_chunk = 1;
  for (int i = 0; i < nthreads; ++i) {
    if (graph->chunkcnt[i] > 1) {
      single_chunk = 0;
      break;
    }
  }
  if (single_chunk) {
    graph->free_adjncy = 1;
    graph->free_adjwgt = 1;

    {
      time_t current_time;
      struct tm * time_info;
      char time_string[9];
      current_time = time(NULL);
      time_info = localtime(&current_time);
      strftime(time_string, 9, "%H:%M:%S", time_info);
      fprintf(stderr, "[%s] single chunk contraction\n", time_string);
    }

    // the chunk structure is backward-compatible when there is only one chunk
    par_contract_graph(ctrl,graph,cnvtxs,(vtx_type const **)gmatch,fcmap);
    // graph->coarser->free_adjncy = 1;
    // graph->coarser->free_adjwgt = 1;

    {
      // manually maintain the chunk info structures
      graph->coarser->chunkcnt[myid] = 1;
      graph->coarser->chunkofst[myid] = vtx_alloc(2);
      graph->coarser->chunkofst[myid][0] = 0;
      graph->coarser->chunkofst[myid][1] = graph->coarser->mynvtxs[myid];
    }
  } else {
    {
      time_t current_time;
      struct tm * time_info;
      char time_string[9];
      current_time = time(NULL);
      time_info = localtime(&current_time);
      strftime(time_string, 9, "%H:%M:%S", time_info);
      fprintf(stderr, "[%s] multi chunk contraction\n", time_string);
    }

    par_contract_chunk_graph(ctrl,graph,cnvtxs,(vtx_type const **)gmatch,fcmap);
  }

  dlthread_free_shmem(gmatch,ctrl->comm);

  dl_free(fcmap);
  dl_free(match);

  // DL_ASSERT(check_graph(graph->coarser),"Invalid graph");

  if (myid == 0) {
    dl_stop_timer(&ctrl->timers.coarsening);
  }

  return graph->coarser;
}




#endif
