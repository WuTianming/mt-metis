/**
 * @file project.c
 * @brief Projection functions. 
 * @author Dominique LaSalle <mtmetis@domnet.org>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-05-20
 */




#ifndef PROJECT_C
#define PROJECT_C




#include "project.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Project a kway partitioning.
 *
 * @param ctrl The control structure containing runtime parameters.
 * @param graph The partitioned graph to project the partition from.
 */
static void S_project_kway(
    ctrl_type * const ctrl,
    graph_type * const graph)
{
  vtx_type i, k, pi, lvtx, nbrid;
  adj_type j, l, istart, iend;
  pid_type mepart, otherpart, na;
  wgt_type tid, ted;
  kwnbrinfo_type * myrinfo;
  adjinfo_type * mynbrs;
  vtx_type * id, * ed;
  pid_type * htable;
  vtx_iset_t * bnd;
  kwnbrinfo_type * nbrinfo;
  kwinfo_type * gkwinfo, * kwinfo;

  tid_type const myid = dlthread_get_id(ctrl->comm);

  pid_type const nparts = ctrl->nparts;
  graph_type * const cgraph = graph->coarser;
  vtx_type * const * const gcmap = graph->cmap;
  size_t   const * const gchunkcnt = graph->chunkcnt;
  vtx_type const * const * const gchunkofst = (vtx_type const **)graph->chunkofst;
  adj_type const * const * const gxadj = (adj_type const **)graph->xadj;
  vtx_type * const * const gadjncy = (vtx_type **)graph->adjncy;
  wgt_type * const * const gadjwgt = (wgt_type **)graph->adjwgt;

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  pid_type ** const gwhere = graph->where;
  pid_type * const where = gwhere[myid];

  vtx_type const mynvtxs = graph->mynvtxs[myid];

  adj_type const * const xadj = gxadj[myid];

  gkwinfo = cgraph->kwinfo;
  dlthread_barrier(ctrl->comm);   // There was a race cond. in the original code
  
  if (myid == 0) {
    wgt_copy(graph->pwgts,cgraph->pwgts,nparts);

    /* migrate old kwinfo */
    graph->kwinfo = gkwinfo;
    cgraph->kwinfo = NULL;
  }
  dlthread_barrier(ctrl->comm);

  kwinfo = gkwinfo+myid;

  /* Compute the required info for refinement */
  kwinfo->nnbrpool = 0;

  /* expanpd nbrinfo */
  dl_free(kwinfo->nbrinfo);
  vtx_iset_free(kwinfo->bnd);
  nbrinfo = kwinfo->nbrinfo = kwnbrinfo_alloc(mynvtxs);
  
  for (i=0;i<mynvtxs;++i) {
    myrinfo = nbrinfo + i;
    myrinfo->nnbrs = 0;
    myrinfo->id = 0;
    myrinfo->ed = 0;
    myrinfo->nbrstart = NULL_ADJ;
  }

  bnd = kwinfo->bnd = vtx_iset_create(0,mynvtxs);

  htable = pid_init_alloc(NULL_PID,nparts);

  // 如果 i 对应的 cnode 拥有跨出去的边，则它属于边界结点，加到 ed 里面去。
  // ed 的比例一开始比较高（90%），逐渐减少，最后减少到50%左右

  // 目标是减少磁盘读写，所以直接不区分 nid 和 ned，完全按照边的存储顺序来。最后处理出所有的 id 和 ed 的值
  // 还要处理出：所有边界结点的 mynbrs 集合。判定 is_bnd()

  int single_chunk = 1;
  for (int i = 0; i < dlthread_get_nthreads(ctrl->comm); ++i) {
    if (graph->chunkcnt[i] > 1) {
      single_chunk = 0;
      break;
    }
  }

  char fin1[1024], fin2[1024];
  sprintf(fin1, "dump_graph%zu_dadjncy_tid%"PF_TID_T".bin", graph->level, myid);
  sprintf(fin2, "dump_graph%zu_dadjwgt_tid%"PF_TID_T".bin", graph->level, myid);
  FILE *adjncy_read = fopen(fin1, "rb");
  FILE *adjwgt_read = fopen(fin2, "rb");
  DL_ASSERT(single_chunk || adjncy_read != NULL, "open adjncy file for read");
  DL_ASSERT(single_chunk || adjwgt_read != NULL, "open adjwgt file for read");

  for (int c=0; c<gchunkcnt[myid]; ++c) {
    vtx_type cstart = gchunkofst[myid][c], cend = gchunkofst[myid][c+1];
    adj_type cadjstart = xadj[cstart], cadjend = xadj[cend];
    if (!single_chunk) {
      fseek(adjncy_read, sizeof(vtx_type) * cadjstart, SEEK_SET);
      fseek(adjwgt_read, sizeof(wgt_type) * cadjstart, SEEK_SET);
      fread(gadjncy[myid], sizeof(vtx_type), cadjend - cadjstart, adjncy_read);
      fread(gadjwgt[myid], sizeof(wgt_type), cadjend - cadjstart, adjwgt_read);
    }
    vtx_type const *const padjncy = gadjncy[myid] - cadjstart;
    wgt_type const *const padjwgt = gadjwgt[myid] - cadjstart;

    for (i=cstart;i<cend;++i) {
      istart = xadj[i];
      iend = xadj[i+1];

      myrinfo = nbrinfo+i;

      DL_ASSERT_EQUALS(myrinfo->nbrstart,NULL_ADJ,"%"PF_ADJ_T);

      /**
       * Caveat NOTE: kwinfo_nbrs can get as big as {EDGES}!!!!!
       * 
       * >>> size(kwinfo_nbrs) = O(min(#part * |V|, |E|))
       * Partitioning into 8 parts and 512 parts are very different!
       */
      na = dl_min(nparts,xadj[i+1]-xadj[i]);
      mynbrs = kwinfo_get_nbrs(kwinfo,i,na);  // list of all neighbouring parts.

      mepart = where[i];
      tid = 0;
      ted = 0;
      for (j=istart; j<iend; ++j) {
        k = padjncy[j];
        if (k < mynvtxs) {
          lvtx = k;
          nbrid = myid;
        } else {
          lvtx = gvtx_to_lvtx(k,graph->dist);
          nbrid = gvtx_to_tid(k,graph->dist);
        }
        otherpart = gwhere[nbrid][lvtx];
        if (mepart == otherpart) {
          tid += padjwgt[j];
        } else {
          ted += padjwgt[j];
          if ((l = htable[otherpart]) == NULL_PID) {
            htable[otherpart] = myrinfo->nnbrs;
            mynbrs[myrinfo->nnbrs].pid = otherpart;
            mynbrs[myrinfo->nnbrs].ed = padjwgt[j];
            ++myrinfo->nnbrs;
          } else {
            mynbrs[l].ed += padjwgt[j];
          }
        }
        DL_ASSERT(myrinfo->nnbrs <= na,"Maxnnbrs = %"PF_PID_T", nnbrs = %" \
            PF_PID_T" for vertex %"PF_TID_T":%"PF_VTX_T"\n",na,myrinfo->nnbrs, \
            myid,i);
      }
      myrinfo->id = tid;
      myrinfo->ed = ted;

      if (ted > 0) {
        if (is_bnd(tid,ted,greedy)) {
          vtx_iset_add(i,bnd);
        }
        for (j=0; j<myrinfo->nnbrs; ++j) {
          htable[mynbrs[j].pid] = NULL_ADJ;
        }
      } else if (tid == 0) {
        /* keep islands on the border */
        vtx_iset_add(i,bnd);
      }
      if (myrinfo->nnbrs == 0) {
        // revoke the neighbor list space just allocated to i
        kwinfo->nnbrpool -= na;
        myrinfo->nbrstart = NULL_ADJ;
      }
    }
  }

  dl_free(htable);

  if (!single_chunk) {
    fclose(adjncy_read);
    fclose(adjwgt_read);
  }

  DL_ASSERT((dlthread_barrier(ctrl->comm),check_kwinfo(kwinfo,graph, \
          (pid_type const **)gwhere)),"Bad info");

  dlthread_barrier(ctrl->comm);
}


/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


void par_project_chunk_graph(
    ctrl_type * const ctrl,
    graph_type * const graph)
{
  vtx_type i, k, lvtx;
  tid_type nbrid;
  pid_type * where;

  tid_type const myid = dlthread_get_id(ctrl->comm);

  vtx_type const * const cmap = graph->cmap[myid];

  vtx_type const mynvtxs = graph->mynvtxs[myid];
  graph_type const * const cgraph = graph->coarser;
  vtx_type const mycnvtxs = cgraph->mynvtxs[myid];
  pid_type const * const * const gcwhere = (pid_type const **)cgraph->where;

  dlthread_barrier(ctrl->comm);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.projection));
  }

  par_graph_alloc_partmemory(ctrl,graph);

  where = graph->where[myid];

  for (i=0;i<mynvtxs;++i) {
    k = cmap[i];
    if (k < mycnvtxs) {
      lvtx = k;
      nbrid = myid;
    } else {
      lvtx = gvtx_to_lvtx(k,cgraph->dist);
      nbrid = gvtx_to_tid(k,cgraph->dist);
    }
    where[i] = gcwhere[nbrid][lvtx];
  }

  /* project refinement information */
  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_ND:
    case MTMETIS_PTYPE_VSEP:
    case MTMETIS_PTYPE_RB:
    case MTMETIS_PTYPE_ESEP:
      dl_error("Unsupported partition type '%d' for chunk-based projection\n",ctrl->ptype);
      break;

    case MTMETIS_PTYPE_KWAY:
      graph->mincut = cgraph->mincut;
      S_project_kway(ctrl,graph);
      break;

    default:
      dl_error("Unknown partition type '%d'\n",ctrl->ptype);
  }

  par_graph_free(graph->coarser);

  if (myid == 0) {
    graph->coarser = NULL;
    dl_stop_timer(&(ctrl->timers.projection));
  }
}




#endif
