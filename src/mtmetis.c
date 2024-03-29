/**
 * @file mtmetis.c
 * @brief Library entry points
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-01
 */




#ifndef MTMETIS_C
#define MTMETIS_C




#include "base.h"
#include "ctrl.h"
#include "graph.h"
#include "partition.h"
#include "order.h"

#include <sys/mman.h>



/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct arg_type {
  ctrl_type * ctrl;
  vtx_type nvtxs;
  vtx_type ncon;
  adj_type const * xadj;
  vtx_type const * adjncy;
  wgt_type const * vwgt;
  wgt_type const * adjwgt;
  int is_mmaped;
  pid_type * where;
  wgt_type * r_obj;
} arg_type;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void S_unify_where(
    pid_type * const * const where, 
    graph_type const * const graph, 
    pid_type * const uwhere)
{
  vtx_type i;

  tid_type const myid = dlthread_get_id(graph->comm);

  for (i=0;i<graph->mynvtxs[myid];++i) {
    uwhere[graph->label[myid][i]] = where[myid][i];
  }
}


static void S_launch_func(
    void * const ptr)
{
  arg_type * arg;
  ctrl_type * ctrl;
  graph_type * graph;
  pid_type * where;
  pid_type ** dwhere;
  tid_type myid, nthreads;

  arg = (arg_type*)ptr;
  ctrl = arg->ctrl;
  where = arg->where;

  myid = dlthread_get_id(ctrl->comm);
  nthreads = dlthread_get_nthreads(ctrl->comm);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.preprocess));
  }

  dwhere = dlthread_get_shmem(sizeof(*dwhere)*nthreads,ctrl->comm);

  /* distribute graph */
  // graph = par_graph_distribute(ctrl->dist,arg->nvtxs,arg->xadj, \
  //     arg->adjncy,arg->vwgt,arg->adjwgt,ctrl->comm);
  graph = par_graph_distribute(ctrl->dist,arg->nvtxs,arg->ncon,ctrl->adjchunksize,arg->xadj, \
      arg->adjncy,arg->is_mmaped,arg->vwgt,arg->adjwgt,ctrl->comm);

  // optimization:
  // free those from wildriver (read in)

  if (myid == 0) {
    if (arg->xadj) {
      if (arg->is_mmaped)
        munmap(arg->xadj, sizeof(adj_type) * (graph->nvtxs+1));
      else
        dl_free(arg->xadj);
    }
    if (arg->adjncy) {
      if (arg->is_mmaped)
        munmap(arg->adjncy, sizeof(vtx_type) * graph->nedges);
      else
        dl_free(arg->adjncy);
    }
    if (arg->vwgt) {
      if (arg->is_mmaped)
        munmap(arg->vwgt, sizeof(wgt_type) * graph->nvtxs * arg->ncon);
      else
        dl_free(arg->vwgt);
    }
    if (arg->adjwgt) {
      dl_free(arg->adjwgt);
    }
  }
  dlthread_barrier(ctrl->comm);

  /* allocate local output vector */
  dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.preprocess));
  }

  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_RB:
      par_partition_rb(ctrl,graph,dwhere);
      break;
    case MTMETIS_PTYPE_KWAY:
      par_partition_kway(ctrl,graph,dwhere);
      break;
    case MTMETIS_PTYPE_ESEP:
      par_partition_edgeseparator(ctrl,graph,dwhere);
      break;
    case MTMETIS_PTYPE_VSEP:
      if (nthreads > 1) {
        par_partition_pre(ctrl,graph);
        /* resize the dwhere allocation to handle the resized graph */
        dwhere[myid] = pid_realloc(dwhere[myid], graph->mynvtxs[myid]);
      }
      par_partition_vertexseparator(ctrl,graph,dwhere);
      break;
    case MTMETIS_PTYPE_ND:
      par_order_nd(ctrl,graph,dwhere);
      break;
    default:
      dl_error("Unknown ptype '%d'",ctrl->ptype);
  }

  if (ctrl->ptype != MTMETIS_PTYPE_ND && \
      ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dlthread_barrier(ctrl->comm);
    if (myid == 0) {
      partition_print_info(ctrl,graph,(pid_type const **)dwhere);
    }
    dlthread_barrier(ctrl->comm);
  }

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.postprocess));
  }

  if (myid == 0 && arg->r_obj) {
    switch (ctrl->ptype) {
      case MTMETIS_PTYPE_ND:
      case MTMETIS_PTYPE_VSEP:
        *(arg->r_obj) = graph->minsep;
        break;
      case MTMETIS_PTYPE_ESEP:
      case MTMETIS_PTYPE_RB:
      case MTMETIS_PTYPE_KWAY:
        *(arg->r_obj) = graph->mincut;
        break;
      default:
        dl_error("Unknown partition type '%d'\n",ctrl->ptype);
    }
  }

  if (where) {
    S_unify_where(dwhere,graph,where);
  }

  dl_free(dwhere[myid]);
  par_graph_free(graph);

  dlthread_free_shmem(dwhere,ctrl->comm);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.postprocess));
  }
}



static char const * S_bool2str(
    int const b)
{
  if (b) {
    return "Yes";
  } else {
    return "No";
  }
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


double * mtmetis_init_options(void)
{
  double * options = double_init_alloc(MTMETIS_VAL_OFF,MTMETIS_NOPTIONS);
  return options;
}


int mtmetis_partition_explicit(
    vtx_type const nvtxs,
    int ncon,
    adj_type const * const xadj,
    vtx_type const * const adjncy,
    int is_mmaped,
    wgt_type const * vwgt,
    wgt_type const * adjwgt,
    double const * const options,
    pid_type * const where,
    wgt_type * const r_objective)
{
  int rv;
  arg_type arg;
  timers_type * timers;
  ctrl_type * ctrl = NULL;
  pid_type ** dwhere = NULL;
  
  if ((rv = ctrl_parse(options,&ctrl)) != MTMETIS_SUCCESS) {
    goto CLEANUP;
  }
  
  ctrl_setup(ctrl,NULL,nvtxs);
  ctrl->maxvwgt = wgt_alloc(ncon);

  timers = &(ctrl->timers);

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dl_print_header("PARAMETERS",'%');
    printf("Number of Threads: %"PF_TID_T" | Verbosity: %s\n",ctrl->nthreads, \
        trans_verbosity_string(ctrl->verbosity));
    printf("Number of Runs: %zu | Random Seed: %u\n",ctrl->nruns, \
        ctrl->seed);
    printf("Number of Partitions: %"PF_PID_T" | Partition Type: %s\n", \
        ctrl->nparts,trans_ptype_string(ctrl->ptype));
    printf("Coarsening Type: %s | Contraction Type: %s\n", \
        trans_ctype_string(ctrl->ctype),trans_contype_string(ctrl->contype));
    printf("Refinement Type: %s | Number of Refinement Passes: %zu\n",
        trans_rtype_string(ctrl->rtype),ctrl->nrefpass);
    printf("Balance: %0.2lf | Distribution: %s\n",ctrl->ubfactor, \
        trans_dtype_string(ctrl->dist));
    printf("Leaf-Matching: %s | Remove Islands: %s\n", \
        S_bool2str(ctrl->leafmatch),S_bool2str(ctrl->removeislands));
    printf("On-Disk: %s\n", \
        S_bool2str(ctrl->ondisk));
    dl_print_footer('%');
  }

  dl_start_timer(&(timers->total));

  arg.ctrl = ctrl;
  arg.nvtxs = nvtxs;
  arg.ncon = ncon;
  arg.xadj = xadj;
  arg.adjncy = adjncy;
  arg.is_mmaped = is_mmaped;
  arg.vwgt = vwgt;
  if (adjwgt && ( \
        ctrl->ptype == MTMETIS_PTYPE_ND || \
        ctrl->ptype == MTMETIS_PTYPE_VSEP)) {
    arg.adjwgt = NULL;
  } else {
    arg.adjwgt = adjwgt;
  }
  arg.where = where;
  arg.r_obj = r_objective;

  dlthread_launch(ctrl->nthreads,&S_launch_func,&arg);

  if (ctrl->time) {
    dl_print_header("MTMETIS TIME",'$');
    printf("Total Time: %.03fs\n",dl_poll_timer(&(timers->total)));
    printf("\tPreprocessing: %.05fs\n",dl_poll_timer(&(timers->preprocess)));
    if (ctrl->ptype == MTMETIS_PTYPE_ND) {
      printf("\tOrdering: %.05fs\n",dl_poll_timer(&(timers->ordering)));
    }
    printf("\tPartitioning: %.05fs\n",dl_poll_timer(&(timers->partitioning)));
    printf("\t\tCoarsening: %.05fs\n",dl_poll_timer(&(timers->coarsening)));
    printf("\t\t\tMatching: %.05fs\n",dl_poll_timer(&(timers->matching)));
    printf("\t\t\tContraction: %.05fs\n", \
        dl_poll_timer(&(timers->contraction)));
    printf("\t\tInitial Partitioning: %.05fs\n",
        dl_poll_timer(&(timers->initpart))); \
    printf("\t\tUncoarsening: %.05fs\n", \
        dl_poll_timer(&(timers->uncoarsening)));
    printf("\t\t\tProjection: %.05fs\n",dl_poll_timer(&(timers->projection)));
    printf("\t\t\tRefinement: %.05fs\n",dl_poll_timer(&(timers->refinement)));
    if (ctrl->ptype == MTMETIS_PTYPE_ND || ctrl->ptype == MTMETIS_PTYPE_RB) {
      printf("\t\tRecursion Overhead: %.05fs\n", \
          dl_poll_timer(&(timers->recursion)));
    }
    if (ctrl->ptype == MTMETIS_PTYPE_ND || ctrl->ptype == MTMETIS_PTYPE_RB || \
        ctrl->metis_serial) {
      printf("\tMetis: %.05fs\n",dl_poll_timer(&(timers->metis)));
    }
    printf("\tPostprocessing: %.05fs\n",dl_poll_timer(&(timers->postprocess)));
    dl_print_footer('$');
  }

  if (ctrl->runstats) {
    dl_print_header("STATISTICS",'&');
    printf("Best Objective: %"PF_WGT_T"\n",wgt_min_value(ctrl->runs, \
          ctrl->nruns));
    printf("Worst Objective: %"PF_WGT_T"\n",wgt_max_value(ctrl->runs, \
          ctrl->nruns));
    printf("Median Objective: %"PF_WGT_T"\n",wgt_median(ctrl->runs, \
          ctrl->nruns));
    printf("Mean Objective - Geo.: %0.2lf - Ari.: %.2lf\n", \
        wgt_geometric_mean(ctrl->runs,ctrl->nruns), \
        wgt_arithmetic_mean(ctrl->runs,ctrl->nruns));
    dl_print_footer('&');
  }


  CLEANUP:

  if (ctrl) {
    ctrl_free(ctrl);
  }
  if (dwhere) {
    r_pid_free(dwhere,ctrl->nthreads);
  }

  return rv;
}



/******************************************************************************
* METIS REPLACEMENTS **********************************************************
******************************************************************************/


int MTMETIS_PartGraphRecursive(
    vtx_type const * const nvtxs,
    vtx_type const * const ncon,
    adj_type const * const xadj,
    vtx_type const * const adjncy,
    wgt_type const * const vwgt,
    vtx_type const * const vsize,
    wgt_type const * const adjwgt,
    pid_type const * const nparts,
    real_type const * const tpwgts,
    real_type const * const ubvec,
    double const * const options,
    wgt_type * const r_edgecut,
    pid_type * const where)
{
  int rv;
  double * modopts;

  modopts = double_duplicate(options,MTMETIS_NOPTIONS);

  modopts[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_RB;
  modopts[MTMETIS_OPTION_NPARTS] = *nparts;

  rv = mtmetis_partition_explicit(*nvtxs,*ncon,xadj,adjncy,0,vwgt,adjwgt,modopts, \
      where,r_edgecut);

  dl_free(modopts);

  return rv;
}


int MTMETIS_PartGraphKway(
    vtx_type const * const nvtxs,
    vtx_type const * const ncon,
    adj_type const * const xadj,
    vtx_type const * const adjncy,
    wgt_type const * const vwgt,
    vtx_type const * const vsize,
    wgt_type const * const adjwgt,
    pid_type const * const nparts,
    real_type const * const tpwgts,
    real_type const * const ubvec,
    double const * const options,
    wgt_type * const r_edgecut,
    pid_type * const where)
{
  int rv;
  double * modopts;

  modopts = double_duplicate(options,MTMETIS_NOPTIONS);

  modopts[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_KWAY;
  modopts[MTMETIS_OPTION_NPARTS] = *nparts;

  rv = mtmetis_partition_explicit(*nvtxs,*ncon,xadj,adjncy,0,vwgt,adjwgt,modopts, \
      where,r_edgecut);

  dl_free(modopts);

  return rv;
}


int MTMETIS_NodeND(
    vtx_type const * const nvtxs,
    adj_type const * const xadj,
    vtx_type const * const adjncy,
    wgt_type const * const vwgt,
    double const * const options,
    pid_type * const perm,
    pid_type * const iperm)
{
  int rv;
  double * modopts;
  vtx_type i;

  if (!iperm) {
    return MTMETIS_ERROR_INVALIDINPUT;
  }

  modopts = double_duplicate(options,MTMETIS_NOPTIONS);
  
  modopts[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_ND;
  modopts[MTMETIS_OPTION_NPARTS] = 3;

  rv = mtmetis_partition_explicit(*nvtxs,1,xadj,adjncy,0,vwgt,NULL,modopts, \
      iperm,NULL);

  /* generate the inverse permutation */
  if (rv == MTMETIS_SUCCESS && perm) {
    for (i = 0; i < *nvtxs; ++i) {
      perm[iperm[i]] = i;
    }
  }

  dl_free(modopts);

  return rv;
}


int MTMETIS_ComputeVertexSeparator(
    vtx_type const * const nvtxs,
    adj_type const * const xadj,
    vtx_type const * const adjncy,
    wgt_type const * const vwgt,
    double const * const options,
    wgt_type * const sepsize,
    pid_type * const where)
{
  int rv;
  double * modopts;

  modopts = double_duplicate(options,MTMETIS_NOPTIONS);
  
  modopts[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_VSEP;
  modopts[MTMETIS_OPTION_NPARTS] = 3;

  rv = mtmetis_partition_explicit(*nvtxs,1,xadj,adjncy,0,vwgt,NULL,modopts, \
      where, sepsize);

  dl_free(modopts);

  return rv;
}




#endif
