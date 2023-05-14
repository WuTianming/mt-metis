/**
 * @file kwayrefine.c
 * @brief KWay refinement routines
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2013-2015, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_KWAYREFINE_C
#define MTMETIS_KWAYREFINE_C




#include "kwayrefine.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct move_type {
  pid_type to;
  vtx_type vtx;
} move_type;


typedef struct update_type {
  pid_type to;
  pid_type from;
  wgt_type ewgt;
  vtx_type nbr;
} update_type;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX vw
#define DLPQ_KEY_T wgt_type
#define DLPQ_VAL_T vtx_type
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLLCB_PREFIX move
#define DLLCB_TYPE_T move_type
#define DLLCB_STATIC 1
#include "dllcb_headers.h"
#undef DLLCB_STATIC
#undef DLLCB_TYPE_T
#undef DLLCB_PREFIX


#define DLLCB_PREFIX update
#define DLLCB_TYPE_T update_type
#define DLLCB_STATIC 1
#include "dllcb_headers.h"
#undef DLLCB_STATIC
#undef DLLCB_TYPE_T
#undef DLLCB_PREFIX


#define DLPQ_PREFIX vt
#define DLPQ_KEY_T wgt_type
#define DLPQ_VAL_T vtx_type
#define DLPQ_USE_HT
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_USE_HT
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLMSET_PREFIX vtx
#define DLMSET_TYPE_T vtx_type
#define DLMSET_STATIC
#include "dlmset_headers.h"
#undef DLMSET_STATIC
#undef DLMSET_TYPE_T
#undef DLMSET_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const MIN_HILL_SIZE = 3;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline int S_right_side(
    int const dir,
    pid_type const to,
    pid_type const from)
{
  if (dir) {
    return to < from;
  } else {
    return to > from;
  }
}


/**
 * @brief Update a vertex incident the one being moved.
 *
 * @param ctrl The control structure.
 * @param k The vertex to update.
 * @param to The partition the vertex is being moved to.
 * @param from The partition the vertex is being moved from.
 * @param ewgt The edge connecting the moved vertex.
 * @param graph The graph.
 * @param bnd The boundary set of vertices.
 * @param queue The priority queue of boundary vertices to move.
 *
 * @return The change in edgecut.
 */
static wgt_type S_update_vertex(
    ctrl_type * ctrl, 
    vtx_type const k, 
    pid_type const to, 
    pid_type const from, 
    wgt_type const ewgt, 
    graph_type * graph, 
    kwinfo_type * const kwinfo,
    vw_pq_t * queue)
{
  vtx_type l;
  wgt_type oed;
  real_type rgain;
  vtx_iset_t * bnd;
  kwnbrinfo_type * myrinfo;
  adjinfo_type * mynbrs;

  tid_type const myid = dlthread_get_id(ctrl->comm);

  pid_type * const where = graph->where[myid];
  pid_type const nparts = ctrl->nparts;
  pid_type const me = where[k];

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  bnd = kwinfo->bnd;

  /* create my workspace */
  myrinfo = kwinfo->nbrinfo+k;

  oed = myrinfo->ed;
  
  mynbrs = kwinfo_get_nbrs(kwinfo,k, \
      dl_min(nparts,graph->xadj[myid][k+1]-graph->xadj[myid][k]));

  if (me == to) {
    myrinfo->id += ewgt;
    myrinfo->ed -= ewgt;
  } else if (me == from) {
    myrinfo->id -= ewgt;
    myrinfo->ed += ewgt;
  }
  /* add it to the boundary if necessary */
  if (!vtx_iset_contains(k,bnd)) {
    if (is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_add(k,bnd);
    }
  } else if (!is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
    vtx_iset_remove(k,bnd);
  }

  /* update nbrs */
  for (l=0;l<myrinfo->nnbrs;++l) {
    if (mynbrs[l].pid == from) {
      if (mynbrs[l].ed == ewgt) {
        mynbrs[l] = mynbrs[--myrinfo->nnbrs];
      } else {
        mynbrs[l].ed -= ewgt;
      }
      break;
    }
  }
  for (l=0;l<myrinfo->nnbrs;++l) {
    if (mynbrs[l].pid == to) {
      mynbrs[l].ed += ewgt;
      break;
    }
  }
  if (to != me && l == myrinfo->nnbrs) {
    mynbrs[myrinfo->nnbrs].ed = ewgt;
    mynbrs[myrinfo->nnbrs++].pid = to;
  }

  if (queue) {
    rgain = ((myrinfo->nnbrs > 0 ? \
          ((real_type)myrinfo->ed) / sqrt(myrinfo->nnbrs) : 0.0) \
          - myrinfo->id);

    if ((me == to || me == from)) {
      if (vw_pq_contains(k,queue)) {
        if ((!greedy && myrinfo->ed > 0) || rgain >= 0) {
          vw_pq_update(rgain,k,queue);
        } else {
          vw_pq_remove(k,queue);
        }
      }
    } // this function guarantees no node will be added into the pq
  }

  return oed - myrinfo->ed;
}


static wgt_type S_move_vertex(      // uses combuffer_add
    ctrl_type * const ctrl,
    graph_type * const graph,
    adj_type const chunkoffset,
    tid_type const myid,
    vtx_type const i,
    pid_type const to,
    kwinfo_type * const kwinfo,
    wgt_type * const pwgts,
    pid_type * const where,
    vw_pq_t * const q,
    update_combuffer_t * const combuffer)
{
  int ncon = graph->ncon;

  vtx_type k;
  adj_type j;
  wgt_type cut, ted, ewgt;
  tid_type nbrid;
  update_type up;
  adjinfo_type * mynbrs;

  pid_type const nparts = ctrl->nparts;
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  adj_type const * const xadj = graph->xadj[myid];
  vtx_type const * const adjncy = graph->adjncy[myid] - chunkoffset;
  wgt_type const * const adjwgt = graph->adjwgt[myid] - chunkoffset;
  wgt_type const * const vwgt = graph->vwgt[myid];

  vtx_iset_t * const bnd = kwinfo->bnd;

  kwnbrinfo_type * const myrinfo = kwinfo->nbrinfo+i;
  pid_type const from = where[i];

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  cut = 0;

  for (int k=0;k<ncon;++k) {
    pwgts[to*ncon+k] += vwgt[i*ncon+k];
    pwgts[from*ncon+k] -= vwgt[i*ncon+k];
  }
  // pwgts[to] += vwgt[i];
  // pwgts[from] -= vwgt[i];
  where[i] = to;

  ted = myrinfo->ed;

  mynbrs = kwinfo_get_nbrs(kwinfo,i, \
      dl_min(nparts,graph->xadj[myid][i+1]-graph->xadj[myid][i]));

  for (k=0;k<myrinfo->nnbrs;++k) {
    if (mynbrs[k].pid == to) {
      break;
    }
  }
  if (k==myrinfo->nnbrs) {
    k = NULL_PID;
  }
  
  /* make the move */
  if (k != NULL_PID) {
    myrinfo->ed += myrinfo->id-mynbrs[k].ed;
    dl_swap(myrinfo->id,mynbrs[k].ed);
  } else if (myrinfo->id > 0) {
    myrinfo->ed += myrinfo->id;
    mynbrs[myrinfo->nnbrs].ed = myrinfo->id;
    k = myrinfo->nnbrs++;
    myrinfo->id = 0;
  }

  /* old minus new */
  cut += ted - myrinfo->ed;

  if (mynbrs[k].ed == 0) {
    mynbrs[k] = mynbrs[--myrinfo->nnbrs];
  } else {
    mynbrs[k].pid = from;
  }
 
  /* see if this vertex should be removed/added from the boundary */
  if (vtx_iset_contains(i,bnd)) {
    if (!is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_remove(i,bnd);
    }
  } else {
    if (is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_add(i,bnd);
    }
  }

  /* update neighbors */
  for(j=xadj[i];j<xadj[i+1];++j) {
    k = adjncy[j];        // NOTE: 每次挪动一个节点的时候，要通知他的相邻节点更新他们的 id 和 ed
    ewgt = adjwgt[j];
    if (k < mynvtxs) {
      /* I own it */
      cut += S_update_vertex(ctrl,k,to,from,ewgt,graph,kwinfo,q);
    } else {
      /* notify my neighbor */
      nbrid = gvtx_to_tid(k,graph->dist);

/*
      if (graph->xadj[nbrid][gvtx_to_lvtx(k,graph->dist)] == graph->xadj[nbrid][gvtx_to_lvtx(k,graph->dist)+1]) {
        dl_error("Error: I am connected to a node who has no outgoing edges. This is impossible.\n");
      }
*/

      up.to = to;
      up.from = from;
      up.ewgt = ewgt;
      up.nbr = gvtx_to_lvtx(k,graph->dist);

      update_combuffer_add(nbrid,up,combuffer);
    }
  }

  return cut;
}


static inline void S_par_sync_pwgts(
    tid_type const myid,
    pid_type const nparts,
    int const ncon,
    wgt_type * const gpwgts,
    wgt_type * const lpwgts,
    dlthread_comm_t const comm)
{
  pid_type p;

  /* turn local pwgts into deltas */
  for (p=0;p<nparts*ncon;++p) {
    lpwgts[p] -= gpwgts[p];
  }

  /* create global deltas */
  wgt_dlthread_sumareduce(lpwgts,nparts*ncon,comm);

  /* set local pwgts to be global pwgts */
  for (p=0;p<nparts*ncon;++p) {
    lpwgts[p] += gpwgts[p];
  }

  dlthread_barrier(comm);

  /* re-sync global pwgts */
  if (myid == 0) {
    for (p=0;p<nparts*ncon;++p) {
      gpwgts[p] = lpwgts[p];
    }
  }

  dlthread_barrier(comm);
}




/******************************************************************************
* REFINEMENT FUNCTIONS ********************************************************
******************************************************************************/


static vtx_type S_par_kwayrefine_GREEDY(
    ctrl_type * const ctrl, 
    graph_type * const graph,
    size_t const niter, 
    kwinfo_type * const kwinfo)
{
  int ncon = graph->ncon;

  vtx_type c, i, k, nmoved;
  adj_type j;
  wgt_type gain, mycut, ewgt;
  wgt_type const * myvwgt;
  pid_type to, from;
  size_t pass;
  real_type rgain;
  wgt_type * lpwgts;
  kwnbrinfo_type * myrinfo;
  adjinfo_type const * mynbrs;
  vtx_iset_t * bnd;
  vw_pq_t * q;
  wgt_type * minwgt, * maxwgt;
  update_type up;
  update_combuffer_t * combuffer;

  tid_type const myid = dlthread_get_id(ctrl->comm);

  /* Link the graph fields */
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  wgt_type const * const vwgt = graph->vwgt[myid];

  pid_type ** const gwhere = graph->where;
  wgt_type * const pwgts = graph->pwgts;
  
  pid_type const nparts = ctrl->nparts;

  kwnbrinfo_type * const nbrinfo = kwinfo->nbrinfo;
  pid_type * const where = gwhere[myid];

  // the balance factors `tpwgts` are the same for all constraints for now, hence the vector length being `nparts` instead of `nparts * ncon`
  real_type const * const tpwgts = ctrl->tpwgts;

  combuffer = update_combuffer_create(graph->mynedges[myid],ctrl->comm);

  lpwgts = wgt_alloc(nparts * ncon);
  // copy pwgts into lpwgts, length = nparts * ncon
  // part_weights, local_part_weights
  wgt_copy(lpwgts,pwgts,nparts * ncon);

  minwgt = wgt_alloc(nparts * ncon);
  maxwgt = wgt_alloc(nparts * ncon);


  bnd = kwinfo->bnd;

  /* setup max/min partition weights */
  for (i=0;i<nparts;++i) {
    for (int t = 0; t < ncon; ++t) {
      maxwgt[i*ncon+t] = ctrl->tpwgts[i]*graph->tvwgt[t]*ctrl->ubfactor;
      minwgt[i*ncon+t] = ctrl->tpwgts[i]*graph->tvwgt[t]*(1.0/ctrl->ubfactor);
    }
  }

  DL_ASSERT(check_kwinfo(kwinfo,graph,(pid_type const **)gwhere),"Bad kwinfo");
  DL_ASSERT(check_kwbnd(kwinfo->bnd,graph,1),"Bad boundary");

  int single_chunk = 1;
  int maxchunkcnt = 1;
  for (int i = 0; i < dlthread_get_nthreads(ctrl->comm); ++i) {
    if (graph->chunkcnt[i] > 1) {
      single_chunk = 0;
      dl_storemax(maxchunkcnt, graph->chunkcnt[i]);
    }
  }

  char fin1[1024], fin2[1024];
  sprintf(fin1, "dump_graph%zu_dadjncy_tid%"PF_TID_T".bin", graph->level, myid);
  sprintf(fin2, "dump_graph%zu_dadjwgt_tid%"PF_TID_T".bin", graph->level, myid);
  FILE *adjncy_read = fopen(fin1, "rb");
  FILE *adjwgt_read = fopen(fin2, "rb");
  DL_ASSERT(single_chunk || adjncy_read != NULL, "open adjncy file for read");
  DL_ASSERT(single_chunk || adjwgt_read != NULL, "open adjwgt file for read");

  q = vw_pq_create(0,mynvtxs); 

  nmoved = 0;
  vtx_type * cperm = vtx_alloc(maxchunkcnt);
  vtx_incset(cperm, 0, 1, maxchunkcnt);

  for (pass=0; pass<niter; pass++) {
    // shuffle chunk order every time
    unsigned seed = ctrl->seed + myid;
    vtx_shuffle_r(cperm, graph->chunkcnt[myid], &seed);

    int no_improvement = 1;
    for (int cp=0; cp < maxchunkcnt; ++cp) {
      int cc = cperm[cp];   // chunk order is randomly shuffled
      int c1 = dl_min(cc,  graph->chunkcnt[myid]),
          c2 = dl_min(cc+1,graph->chunkcnt[myid]);

      vtx_type cstart = graph->chunkofst[myid][c1],
               cend   = graph->chunkofst[myid][c2];
      adj_type cadjstart = graph->xadj[myid][cstart],
               cadjend   = graph->xadj[myid][cend];

      // read the edge files from the disk
      if (!single_chunk) {
        fseek(adjncy_read, sizeof(vtx_type) * cadjstart, SEEK_SET);
        fseek(adjwgt_read, sizeof(wgt_type) * cadjstart, SEEK_SET);
        fread(graph->adjncy[myid], sizeof(vtx_type), cadjend - cadjstart, adjncy_read);
        fread(graph->adjwgt[myid], sizeof(wgt_type), cadjend - cadjstart, adjwgt_read);
      }

      mycut = 0;
      for (c=0;c<2;++c) {
        dlthread_barrier(ctrl->comm);

        /* fill up my queue with my vertices */
        /* in the chunk-based version, only vertices within current chunk are inserted */
        vw_pq_clear(q);
        if (c1 < c2) for (i=0;i<bnd->size;++i) {
          k = vtx_iset_get(i,bnd);
          if (k < cstart || k >= cend) continue;    // only insert current chunk nodes
          DL_ASSERT(k < mynvtxs,"Invalid border vertex %"PF_VTX_T,k);
          if (nbrinfo[k].nnbrs > 0) {
            /* only insert vertices with external neighbors */
            rgain = (1.0*nbrinfo[k].ed/sqrt(nbrinfo[k].nnbrs)) - nbrinfo[k].id;
            vw_pq_push(rgain,k,q);
          }
        }
        /* make moves */

        do {
          /* perform updates */
          while (update_combuffer_next(&up,combuffer)) {
            k = up.nbr;
            ewgt = up.ewgt;
            to = up.to;
            from = up.from;

            mycut += S_update_vertex(ctrl,k,to,from,ewgt,graph,kwinfo,q);
          }

          /* move a vertice */
          if (q->size > 0) {
            i = vw_pq_pop(q);

            myrinfo = kwinfo->nbrinfo+i;
            mynbrs = kwinfo_get_nbrs_ro(kwinfo,i, \
                dl_min(nparts,graph->xadj[myid][i+1]-graph->xadj[myid][i]));

            from = where[i];
            myvwgt = vwgt + i * ncon;

            if (myrinfo->id > 0) {
              int give_up = 0;
              for (int t=0; t<ncon; ++t) {
                if (lpwgts[from*ncon+t]-myvwgt[t] < minwgt[from*ncon+t]) {
                  give_up = 1;
                  break;
                }
              }
              if (give_up) {
                continue;
              }
            }

            /* find the first eligible partition */
            for (k=0;k<myrinfo->nnbrs; ++k) {
              to = mynbrs[k].pid;
              if (!S_right_side(c,to,from)) {
                continue;
              }
              
              int overweight = 0;
              for (int t=0; t<ncon; ++t) {
                if (lpwgts[to*ncon+t]+myvwgt[t] > maxwgt[to*ncon+t]) {
                  overweight = 1;
                  break;
                }
              }

              if (!overweight) {     // keep nicely balanced after move
                if (mynbrs[k].ed >= myrinfo->id) {    // criterion for improvement
                  break;
                }
              }
            }
            if (k == myrinfo->nnbrs) {
              /* if there aren't any eligable partitions, abort */
              // give_up++;
              continue;
            }

            /* see if there is a better one from the eligible one */
            for (j=k+1; j<myrinfo->nnbrs; ++j) {
              to = mynbrs[j].pid;
              if (!S_right_side(c,to,from)) {
                continue;
              }
              if (mynbrs[j].ed >= mynbrs[k].ed) {
                gain = mynbrs[j].ed-myrinfo->id; 
                DL_ASSERT(gain >= 0, "Invalid gain of %"PF_WGT_T,gain);

                int overweight = 0;
                for (int t=0; t<ncon; ++t) {
                  if (lpwgts[to*ncon+t]+myvwgt[t] > maxwgt[to*ncon+t]) {
                    overweight = 1;
                    break;
                  }
                }

                if ((gain > 0 && !overweight) \
                    || (mynbrs[j].ed == mynbrs[k].ed && \
                      tpwgts[mynbrs[k].pid]*lpwgts[to] < \
                      tpwgts[to]*lpwgts[mynbrs[k].pid])) {
                  k = j;
                }
              }
            }
            to = mynbrs[k].pid;

            if (mynbrs[k].ed >= myrinfo->id) { 
              gain = mynbrs[k].ed-myrinfo->id;

              int overweight1 = 0, overweight2 = 0;
              for (int t=0; t<ncon; ++t) {
                if (lpwgts[from*ncon+t] >= maxwgt[from*ncon+t]) {
                  overweight1 = 1;
                  break;
                }
              }
              for (int t=0; t<ncon; ++t) {
                if (tpwgts[to]*lpwgts[from*ncon+t] > \
                    tpwgts[from]*(lpwgts[to*ncon+t]+myvwgt[t])) {
                  overweight2 = 1;
                  break;
                }
              }

              if (!(gain > 0 || (gain == 0 \
                        && (overweight1 || overweight2)))) {
                continue;
              }
            }

            int overweight = 0;
            int underweight = 0;

            for (int t=0; t<ncon; ++t) {
              if (lpwgts[to*ncon+t]+myvwgt[t] > maxwgt[to*ncon+t]) {
                overweight = 1;
                break;
              }
            }

            for (int t=0; t<ncon; ++t) {
              if (lpwgts[from*ncon+t]-myvwgt[t] < minwgt[from*ncon+t]) {
                underweight = 1;
                break;
              }
            }

            if (overweight || underweight) {
              /* whatever you do, don't push the red button */
              continue;
            }

            /* make the move ***************************************************/
            ++nmoved;

            mycut += S_move_vertex(ctrl,graph,cadjstart,myid,i,to,kwinfo,lpwgts, \
                where,q,combuffer);
          } 
          // OPTIMIZE: should I add some delay here to save CPU?
        } while ((q->size > 0 && vw_pq_top(q) >= 0) || \
            !update_combuffer_finish(combuffer));

        DL_ASSERT_EQUALS(update_combuffer_next(NULL,combuffer),0,"%d");

        update_combuffer_clear(combuffer);

        /* update my partition weights */
        S_par_sync_pwgts(myid,nparts,ncon,pwgts,lpwgts,ctrl->comm);

      } /* end directions */

      mycut = wgt_dlthread_sumreduce(mycut,ctrl->comm);

      par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH, \
          "Refinement pass %zu (%d/%d): %"PF_WGT_T" improvement\n",pass,cp+1,maxchunkcnt,mycut);

      if (mycut != 0) {
        no_improvement = 0;
      }

      if (myid == 0) {
        graph->mincut -= (mycut/2);
      }
    } /* end chunks */

    if (no_improvement) {
      break;
    }
  } /* end passes */

  dl_free(cperm);

  if (!single_chunk) {
    fclose(adjncy_read);
    fclose(adjwgt_read);
  }

  nmoved = vtx_dlthread_sumreduce(nmoved,ctrl->comm);

  vw_pq_free(q);

  dl_free(minwgt);
  dl_free(maxwgt);
  dl_free(lpwgts);

  DL_ASSERT(check_kwinfo(kwinfo,graph,(pid_type const **)gwhere),"Bad kwinfo");
  DL_ASSERT(check_kwbnd(kwinfo->bnd,graph,1),"Bad boundary");
  DL_ASSERT(graph->mincut >= 0,"Invalid mincut of %"PF_WGT_T, \
      graph->mincut);
  DL_ASSERT_EQUALS(graph_cut(graph,(pid_type const**)gwhere), \
      graph->mincut,"%"PF_WGT_T);

  /* implicit barrier */
  update_combuffer_free(combuffer);

  return nmoved;
}



/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


vtx_type par_kwayrefine_chunk(
    ctrl_type * const ctrl,
    graph_type * const graph,
    kwinfo_type * const kwinfo)
{
  vtx_type nmoves; 

  nmoves = 0;

  switch (ctrl->rtype) {
    case MTMETIS_RTYPE_GREEDY:
      nmoves = S_par_kwayrefine_GREEDY(ctrl,graph,ctrl->nrefpass,kwinfo);
      break;
    case MTMETIS_RTYPE_HS:
      dl_error("Unsupported refinement type '%d' for chunk-based K-Way partitions.", ctrl->rtype);
      break;
    case MTMETIS_RTYPE_FM:
      /* use KPM version of FM */
      dl_error("Unsupported refinement type '%d' for chunk-based K-Way partitions.", ctrl->rtype);
      break;
    default:
      dl_error("Unsupported refinement type '%d' for K-Way partitions.", ctrl->rtype);
  }

  par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"%zu) [%"PF_VTX_T" %" \
      PF_ADJ_T"] {%"PF_WGT_T" %"PF_VTX_T"}\n",graph->level,graph->nvtxs, \
      graph->nedges,graph->mincut,nmoves);

  return nmoves;
}




#endif
