/**
 * @file refine.c
 * @brief Refinement functions
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-05-20
 */




#ifndef MTMETIS_REIFNE_C
#define MTMETIS_REFINE_C




#include "refine.h"
#include "kwayrefine.h"
#include "check.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Compute the parameters for performing kway partitioning.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 */
static void S_partparams_kway(
    ctrl_type * const ctrl,
    graph_type * const graph)
{
  vtx_type other,me,i,k,lvtx,nbrid,na;
  adj_type j, l;
  wgt_type mincut;
  kwnbrinfo_type * nbrinfo; 
  kwnbrinfo_type * myrinfo;
  adjinfo_type * mynbrs;
  vtx_iset_t * bnd;
  kwinfo_type * kwinfo;
  wgt_type * mypwgts;
  wgt_type ** gpwgts;

  tid_type const myid = dlthread_get_id(ctrl->comm);
  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

  pid_type const nparts = ctrl->nparts;
  int const ncon = graph->ncon;

  adj_type const * const * const gxadj = (adj_type const **)graph->xadj;
  vtx_type const * const * const gadjncy = (vtx_type const **)graph->adjncy;
  wgt_type const * const * const gvwgt = (wgt_type const **)graph->vwgt;
  wgt_type const * const * const gadjwgt = (wgt_type const **)graph->adjwgt;

  pid_type const * const * const gwhere  = (pid_type const **)graph->where;

  wgt_type * const pwgts = graph->pwgts;

  vtx_type const mynvtxs = graph->mynvtxs[myid];

  adj_type const * const xadj = gxadj[myid];
  vtx_type const * const adjncy = gadjncy[myid];
  wgt_type const * const vwgt = gvwgt[myid];
  wgt_type const * const adjwgt = gadjwgt[myid];
  pid_type const * const where = gwhere[myid];

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  DL_ASSERT(graph->pwgts != NULL,"Non-allocated pwgts");
  DL_ASSERT(graph->where != NULL,"Non-allocated where");

  par_kwinfo_create(ctrl,graph);

  gpwgts = dlthread_get_shmem(sizeof(wgt_type*)*nthreads,ctrl->comm);

  mypwgts = gpwgts[myid] = wgt_init_alloc(0,nparts*ncon);

  /* reset the size of neighbor infor array */
  kwinfo = graph->kwinfo + myid; 
  nbrinfo = kwinfo->nbrinfo;
  kwinfo->nnbrpool = 0;
  bnd = kwinfo->bnd;

  mincut = 0;
  vtx_iset_clear(bnd);

  /* calculate partition weights */
  for (i=0;i<mynvtxs;++i) {
    for (int j=0; j<ncon; ++j)
      mypwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
    // mypwgts[where[i]] += vwgt[i];
  }
  dlthread_barrier(ctrl->comm);

  if (myid == 0) {
    /* someday I need to parallelize this via dlthreads */
    for (i=0; i<nparts;++i) {
      for (k=0; k<ncon; ++k) {
        // printf("refine_chunk.c line %d: ", __LINE__);
        pwgts[i*ncon+k] = 0;
        for (j=0; j<nthreads;++j) {
          pwgts[i*ncon+k] += gpwgts[j][i*ncon+k];
          // printf("%"PF_WGT_T" ", gpwgts[j][i*ncon+k]);
        }
        // printf("\n");
      }
    }
  }

  /* calculate nbrinfo for vertices */
  for (i=0; i<mynvtxs; ++i) {
    myrinfo = nbrinfo+i;
    myrinfo->nnbrs = 0;
    myrinfo->id = 0;
    myrinfo->ed = 0;
    myrinfo->nbrstart = NULL_ADJ;

    na = dl_min(nparts,xadj[i+1]-xadj[i]);
    mynbrs = kwinfo_get_nbrs(kwinfo,i,na);

    me = where[i];

    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        nbrid = gvtx_to_tid(k,graph->dist);
        lvtx = gvtx_to_lvtx(k,graph->dist); 
      }
      if (me == gwhere[nbrid][lvtx]) {    // if part[u] == part[v]
        myrinfo->id += adjwgt[j];
      } else {
        myrinfo->ed += adjwgt[j];
      }
    }

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      mincut += myrinfo->ed;
      for (j=xadj[i]; j<xadj[i+1]; ++j) {
        k = adjncy[j];
        if (k < mynvtxs) {
          lvtx = k;
          nbrid = myid;
        } else {
          nbrid = gvtx_to_tid(k,graph->dist);
          lvtx = gvtx_to_lvtx(k,graph->dist); 
        }
        other = gwhere[nbrid][lvtx];
        if (me != other) {    // part[u] != part[v]
          for (l=0; l<myrinfo->nnbrs; l++) {
            if (mynbrs[l].pid == other) {   // why not just make `mynbrs` fixed size
              mynbrs[l].ed += adjwgt[j];
              break;
            }
          }
          if (l == myrinfo->nnbrs) {
            mynbrs[l].pid = other;
            mynbrs[l].ed  = adjwgt[j];
            myrinfo->nnbrs++;
          }
        }
      }

      /* Only ed-id>=0 nodes are considered to be in the boundary */
      if (is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
        vtx_iset_add(i,bnd);
      }
      DL_ASSERT(myrinfo->nnbrs > 0,"No neighbors.");
    } else if (myrinfo->id == 0) {
      vtx_iset_add(i,bnd);
    } else {
      myrinfo->nbrstart = NULL_ADJ;
      kwinfo->nnbrpool -= na;
      DL_ASSERT_EQUALS(myrinfo->nnbrs,0,"%"PF_ADJ_T);
    }
  } 

  mincut = wgt_dlthread_sumreduce(mincut,ctrl->comm);

  dl_free(gpwgts[myid]);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    graph->mincut = mincut/2;
  }

  dlthread_free_shmem(gpwgts,ctrl->comm);

  /* the checks */
  DL_ASSERT(check_kwinfo(kwinfo,graph,(pid_type const **)gwhere),"Bad info");
}


/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


vtx_type par_refine_chunk_graph(
    ctrl_type * const ctrl,
    graph_type * const graph)
{
  vtx_type nmoves; 

  tid_type const myid = dlthread_get_id(ctrl->comm);

  nmoves = 0;

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.refinement));
  }

  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_ND:
    case MTMETIS_PTYPE_VSEP:
      dl_error("Unsupported partition type '%d' for chunk-based refinement\n",ctrl->ptype);
      break;
    case MTMETIS_PTYPE_RB:
    case MTMETIS_PTYPE_ESEP:
      dl_error("Unsupported partition type '%d' for chunk-based refinement\n",ctrl->ptype);
      break;
    case MTMETIS_PTYPE_KWAY:
      if (graph->kwinfo == NULL) {
        // NOTE: kwinfo is empty <==> the graph is the initial partition.
        // hence we can spare the effort to make it support multiple chunks

        // if (myid == 0)
        //   for (int i = 0; i < ctrl->nthreads; ++i)
        //     if (graph->chunkcnt[i] > 1)
        //       dl_error("The initial partition is not single-chunk -- aborting.\n");
        dlthread_barrier(ctrl->comm);

        S_partparams_kway(ctrl,graph);
      }
      nmoves = par_kwayrefine_chunk(ctrl,graph,graph->kwinfo+myid);
      break;
    default:
      dl_error("Unknown partition type '%d'\n",ctrl->ptype);
  }

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.refinement));
  }

  return nmoves;
}




#endif
