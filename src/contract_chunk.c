/**
 * @file contract_chunk.c
 * @brief Functions for performing contraction with disk-based adjlist chunks
 * @author wtm
 * @version 1
 */




#ifndef MTMETIS_CONTACT_CHUNK_C
#define MTMETIS_CONTACT_CHUNK_C




#include "contract.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct edge_type {
  vtx_type dst;
  wgt_type wgt;
} edge_type;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


#define DEF_MASK_SIZE (0x2000)
static uint32_t const MASK_SIZE = DEF_MASK_SIZE;
static uint32_t const MASK = DEF_MASK_SIZE-1;
static vtx_type const MASK_MAX_DEG = DEF_MASK_SIZE >> 3;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLSORT_PREFIX edge
#define DLSORT_TYPE_T edge_type
#define DLSORT_COMPARE(a,b) ((a).dst < (b).dst)
#define DLSORT_STATIC
#include "dlsort_headers.h"
#undef DLSORT_STATIC
#undef DLSORT_COMPARE
#undef DLSORT_TYPE_T
#undef DLSORT_PREFIX






/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Reverse the bits of a key for a given mask. Use for hashing into a
 * secondary hash table to prevent collisions.
 *
 * @param n The key to reverse.
 * @param mask The mask of the bits to be reversed.
 *
 * @return The reversed key.
 */
static inline vtx_type S_reverse(
    vtx_type const n, 
    vtx_type const mask)
{
  vtx_type r = vtx_reversebits(n);
  int const mb = vtx_downlog2(mask);
  int const vs = sizeof(vtx_type)*8;
  if (vs >= 2*mb) {
    return (r >> (vs - (2*mb))) & mask;
  } else {
    return r >> (vs - mb);
  }
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Adjust the coarse vertex map given new graph distribution paramters. 
 *
 * @param cmap The coarse vertex map.
 * @param mynvtxs The number of vertices this thread owns.
 * @param olddist The old graph distribution parameters.
 * @param newdist The new graph distribution parameters.
 */
static void S_adjust_cmap(
    vtx_type * const cmap, 
    vtx_type const mynvtxs, 
    graphdist_type const olddist, 
    graphdist_type const newdist)
{
  vtx_type i,k;
  tid_type o;

  for (i=0;i<mynvtxs;++i) {
    if (cmap[i] >= olddist.offset) { /* remote vertex */
      k = gvtx_to_lvtx(cmap[i],olddist);
      o = gvtx_to_tid(cmap[i],olddist);
      cmap[i] = lvtx_to_gvtx(k,o,newdist);
    }
  }
}


/**
 * @brief Perform contraction using a dense vector.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param mycnvtxs The number of coarse vertices owned by this thread.
 * @param gmatch The global match array.
 * @param fcmap The first fine vertex for each coarse vertex.
 */
static void S_par_contract_DENSE(
    ctrl_type * const ctrl, 
    graph_type * const graph, 
    vtx_type const mycnvtxs, 
    vtx_type const * const * const gmatch, 
    vtx_type const * const fcmap)
{
  adj_type cnedges, cnedges_upperlimit, l, maxdeg, j, i;
  tid_type o, t;
  vtx_type v, c, cg, k;
  wgt_type ewgt;
  adj_type * table;
  graph_type * cgraph;
  graphdist_type cdist;

  tid_type const myid = dlthread_get_id(ctrl->comm);
  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  adj_type const * const * const gxadj = (adj_type const * const *)graph->xadj;
  vtx_type const * const * const gadjncy = (vtx_type const * const *)graph->adjncy;
  wgt_type const * const * const gvwgt = (wgt_type const * const *)graph->vwgt;
  wgt_type const * const * const gadjwgt = (wgt_type const * const *)graph->adjwgt;

  vtx_type const * const * const gcmap = (vtx_type const **)graph->cmap;

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  cgraph = par_graph_setup_coarse(graph,mycnvtxs);

  cdist = cgraph->dist;

  S_adjust_cmap(graph->cmap[myid],mynvtxs,graph->dist,cdist);

  /* count possible edges (upper limit) */
  cnedges_upperlimit = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fcmap[c];   // v is the first vertex that makes up coarse c
    o = myid;
    l = 0;
    do {
      l += gxadj[o][v+1] - gxadj[o][v];

      // traverse like a linked list (although the match list is either len=1 or len=2)
      // to accommodate cases of self-matching
      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fcmap[c]));

    // the maximum degree possible for a coarse vtx (Caveat: can this get big?)
    dl_storemax(maxdeg,l);

    // an upper limit for total coarse edges
    cnedges_upperlimit += l;
  }

  adj_type * const mycxadj = cgraph->xadj[myid];
  // over-allocating...
  // TODO: apply the new chunk-based memory model to this function
  vtx_type * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(cnedges_upperlimit);
  wgt_type * const mycvwgt = cgraph->vwgt[myid];
  wgt_type * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(cnedges_upperlimit);

  table = NULL;

  table = adj_init_alloc(NULL_ADJ,graph->gnvtxs);

  cnedges = 0;    // coarse edge count is to be calculated accurately later
  mycxadj[0] = 0;

  dlthread_barrier(ctrl->comm);

  /* set max degree for the coarse graph */
  for (c=0;c<mycnvtxs;++c) {
    cg = lvtx_to_gvtx(c,myid,cdist);
    /* initialize the coarse vertex */
    mycvwgt[c] = 0;     // the v weight is the sum of original vertices

    v = fcmap[c];
    o = myid;
    DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
        graph->dist),"%"PF_TID_T);
    do {    // traverse all original vertices (o in v) for this coarse vertex c
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],cdist),"%"PF_VTX_T);

      /* transfer over vertex stuff from v and u */
      mycvwgt[c] += gvwgt[o][v];

      // this kills locality when accessing the fine graph:
      //   (o, v) = canonical(gmatch[o][v]);

      /* for all edges emanating from (v in o): */
      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = gadjncy[o][j];    // TODO slice
        if (k < graph->mynvtxs[o]) {
          t = o;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        } // edge: (v in o) -> (k in t)
        k = gcmap[t][k];
        if (gvtx_to_tid(k,cdist) == myid) {
          k = gvtx_to_lvtx(k,cdist);
        } // k is now the canonical expression of the coarse vertex that j leads to

        if (k == c || k == cg) {
          /* internal edge */
        } else {
          /* external edge */
          i = table[k];
          ewgt = graph->uniformadjwgt ? 1 : gadjwgt[o][j];  // TODO slice
          if (i == NULL_ADJ) {
            /* new edge */
            mycadjncy[cnedges] = k;     // TODO slice
            mycadjwgt[cnedges] = ewgt;  // TODO slice
            table[k] = cnedges++;       // put a new edge in the dense vector
          } else {
            /* duplicate edge */
            mycadjwgt[i] += ewgt;       // TODO slice
          }
        }
      }

      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(myid == o && v == fcmap[c]));

    /* clear the table */
    for (j = cnedges;j > mycxadj[c];) {
      --j;
      k = mycadjncy[j];     // an std::unordered_map would be perfect here...
      table[k] = NULL_ADJ;
    }

    mycxadj[c+1] = cnedges;
  }

  dl_free(table);

  cgraph->mynedges[myid] = cnedges;

  // TODO: why not???
  // but that's OK anyway, because in the ondisk implementation
  //   the space will be freed soon.
  //graph_readjust_memory(cgraph,adjsize);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  par_graph_setup_twgts(cgraph);
  
  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  // Note: the check.c utilities will be unavailable for quite some time
  // because they require a total rewrite
  // DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}


/**
 * @brief Perform contraction using a single hash table, performing a linear
 * scan on the edge list in the case of a collsion.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param mycnvtxs The number of coarse vertices owned by this thread.
 * @param gmatch The global match array.
 * @param fcmap The first fine vertex for each coarse vertex.
 */
static void S_par_contract_CLS_quadratic(
    ctrl_type * const ctrl, 
    graph_type * const graph, 
    vtx_type const mycnvtxs, 
    vtx_type const * const * const gmatch, 
    vtx_type const * const fcmap)
{
  adj_type cnedges, l, maxdeg, j, i, jj, start;
  tid_type o, t;
  vtx_type v, c, cg, k;
  wgt_type ewgt;
  graph_type * cgraph;
  graphdist_type cdist;
  offset_type * htable;

  tid_type const myid = dlthread_get_id(ctrl->comm);
  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  size_t   const * const gchunkcnt = graph->chunkcnt;
  vtx_type const * const * const gchunkofst = (vtx_type const **)graph->chunkofst;
  adj_type const * const * const gxadj = (adj_type const * const *)graph->xadj;
  vtx_type const * const * const gadjncy = (vtx_type const * const *)graph->adjncy;
  wgt_type const * const * const gvwgt = (wgt_type const * const *)graph->vwgt;
  wgt_type const * const * const gadjwgt = (wgt_type const * const *)graph->adjwgt;

  vtx_type const * const * const gcmap = (vtx_type const **)graph->cmap;

  /* count possible edges */
  cnedges = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fcmap[c];
    o = myid;
    l = 0;
    do {
      l += gxadj[o][v+1] - gxadj[o][v];
      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fcmap[c]));
    dl_storemax(maxdeg,l);
    cnedges += l;
  }

  size_t maxchunkcnt = graph->chunkcnt[myid];
  for (tid_type t = 0; t < dlthread_get_nthreads(graph->comm); ++t) {
    dl_storemax(maxchunkcnt, graph->chunkcnt[t]);
  }

/*
  // I have to implement chunks in this too
  // or just disable this branch
/*
  if (maxdeg > MASK_MAX_DEG) {
    S_par_contract_DENSE(ctrl,graph,mycnvtxs,gmatch,fcmap);
    return;
  }
*/

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  cgraph = par_graph_setup_coarse(graph,mycnvtxs);

  cdist = cgraph->dist;

  S_adjust_cmap(graph->cmap[myid],mynvtxs,graph->dist,cdist);

  adj_type * const mycxadj = cgraph->xadj[myid];
  wgt_type * const mycvwgt = cgraph->vwgt[myid];

  size_t adjncy_chunksize = /* take min */
      ctrl->adjchunksize > gxadj[myid][gchunkofst[myid][1]] * 2
          ? gxadj[myid][gchunkofst[myid][1]] * 2
          : ctrl->adjchunksize * 1.1;

  vtx_type * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(adjncy_chunksize);
  wgt_type * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(adjncy_chunksize);

  size_t   * cchunkcnt = cgraph->chunkcnt; cchunkcnt[myid] = 1;
  size_t   * pchunkcnt = cchunkcnt + myid;
  vtx_type * mycchunkofst = cgraph->chunkofst[myid] = vtx_alloc(maxchunkcnt * 1.5 + 5);
  mycchunkofst[0] = 0;

  htable = offset_init_alloc(NULL_OFFSET,MASK_SIZE);

  cnedges = 0;
  mycxadj[0] = 0;

  dlthread_barrier(ctrl->comm);

  char fout1[1024], fout2[1024];
  sprintf(fout1, "dump_graph%zu_dadjncy_tid%"PF_TID_T".bin", cgraph->level, myid);
  sprintf(fout2, "dump_graph%zu_dadjwgt_tid%"PF_TID_T".bin", cgraph->level, myid);
  FILE *cadjncy_dump = fopen(fout1, "wb");
  FILE *cadjwgt_dump = fopen(fout2, "wb");
  DL_ASSERT(cadjncy_dump != NULL, "open cadjncy file");
  DL_ASSERT(cadjwgt_dump != NULL, "open cadjwgt file");

  char fin1[1024], fin2[1024];
  sprintf(fin1, "dump_graph%zu_dadjncy_tid%"PF_TID_T".bin", graph->level, myid);
  sprintf(fin2, "dump_graph%zu_dadjwgt_tid%"PF_TID_T".bin", graph->level, myid);
  FILE *dadjncy_read = fopen(fin1, "rb");
  FILE *dadjwgt_read = fopen(fin2, "rb");
  DL_ASSERT(dadjncy_read != NULL, "open dadjncy file");
  DL_ASSERT(dadjwgt_read != NULL, "open dadjwgt file");

  /**
   * the quadratic contract scheme works by doing the following loop:
   *
   * for cc1 in range(0, gchunkcnt[myid]):
   *   read_adjncy(cc1 into buf1[myid]);
   *   for cc2 in range(cc1+1, max(gchunkcnt)):
   *     read_adjncy(cc2 into buf2[myid]);
   *     now buf2[] has all adjncy info for chunk number = cc2!
   *     forall u in myvtxs && u in cc1 && u == fcmap[gcmap[u]]:
   *       if match[u] in cc2:
   *         copy adjncy[{xadj[u]}] from shared memory into a buf
   *         ## GIVEN cleanup_match puts cvtxs with the same cc2 together,
   *         do contract for all vertices in cc1 with cc2
   *     barrier();
   *   barrier();
  */

  vtx_type * const local_adjncy = vtx_alloc(adjncy_chunksize);
  wgt_type * const local_adjwgt = wgt_alloc(adjncy_chunksize);

  // vtx_type ** dlocal_adjncy = dlthread_get_shmem(sizeof(vtx_type*)*nthreads, graph->comm);

  c = 0;
  twgt_type adjwgt_sum = 0;

  for (int cc1=0; cc1<maxchunkcnt; ++cc1) {
    dlthread_barrier(graph->comm);    // barrier every time before file read
    if (cc1 >= gchunkcnt[myid]) {
      for (int cc2=cc1; cc2<maxchunkcnt; ++cc2) {
        dlthread_barrier(graph->comm);    // barrier every time before file read
      }
      continue;
    }
    vtx_type cc1start = gchunkofst[myid][cc1], cc1end = gchunkofst[myid][cc1+1];
    adj_type cc1adjstart = gxadj[myid][cc1start], cc1adjend = gxadj[myid][cc1end];
    fseek(dadjncy_read, sizeof(vtx_type) * cc1adjstart, SEEK_SET);
    fseek(dadjwgt_read, sizeof(wgt_type) * cc1adjstart, SEEK_SET);
    fread(local_adjncy, sizeof(vtx_type), cc1adjend - cc1adjstart, dadjncy_read);   // 读的是实际长度
    fread(local_adjwgt, sizeof(wgt_type), cc1adjend - cc1adjstart, dadjwgt_read);

// #define PRINT_LOOP_VAR

#ifdef PRINT_LOOP_VAR
      {
        fprintf(stderr, "#%"PF_TID_T": this chunk is [%"PF_VTX_T", %"PF_VTX_T")\n", myid, cc1start, cc1end);
      }
#endif

    for (int cc2=cc1; cc2<maxchunkcnt; ++cc2) {
      dlthread_barrier(graph->comm);    // barrier every time before file read
#ifdef PRINT_LOOP_VAR
      if (myid == 0) {
        fprintf(stderr, "chunk vector (%d,%d)\n", cc1, cc2);
      }
#endif
      vtx_type cc2start, cc2end;
      adj_type cc2adjstart, cc2adjend;
      if (cc2 < gchunkcnt[myid]) {
        cc2start = gchunkofst[myid][cc2], cc2end = gchunkofst[myid][cc2+1];
        cc2adjstart = gxadj[myid][cc2start], cc2adjend = gxadj[myid][cc2end];
        fseek(dadjncy_read, sizeof(vtx_type) * cc2adjstart, SEEK_SET);
        fseek(dadjwgt_read, sizeof(wgt_type) * cc2adjstart, SEEK_SET);
        fread(gadjncy[myid], sizeof(vtx_type), cc2adjend - cc2adjstart, dadjncy_read);
        fread(gadjwgt[myid], sizeof(wgt_type), cc2adjend - cc2adjstart, dadjwgt_read);  // this space is shared
      }

      // coarse vertices have their info gathered in order
      // start where we left off at coarse idx = c
      for (;c<mycnvtxs;++c) {
        cg = lvtx_to_gvtx(c,myid,cdist);
        /* initialize the coarse vertex */
        mycvwgt[c] = 0;

        v = fcmap[c];       // this is guaranteed to be thread local!
        if (v >= cc1end) {  // I'm done with the cc1 chunk
          break;
        }
        DL_ASSERT(v >= cc1start, "CHECK FAILED: #%"PF_TID_T" regression (fcmap[%"PF_VTX_T"/%"PF_VTX_T"]=%"PF_VTX_T") at line %d in contract_chunk\n", myid, c, mycnvtxs, v, __LINE__);
        o = myid;
        DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],cdist),"%"PF_VTX_T);

        k = gmatch[o][v]; t = myid;
        if (k >= graph->mynvtxs[o]) {
          t = gvtx_to_tid(k, graph->dist);
          k = gvtx_to_lvtx(k, graph->dist);
        } // now (o,v) is paired to (t,k)

        if (k >= gchunkofst[t][cc2+1]) {  // done with cc2 chunk
          break;
        }
        DL_ASSERT(k >= gchunkofst[t][cc2], "CHECK FAILED: regression at line %d in contract_chunk\n", __LINE__);
        DL_ASSERT_EQUALS(gcmap[o][v],gcmap[t][k],"%"PF_VTX_T);
        DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[t][k],cdist),"%"PF_VTX_T);

        DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist), graph->dist),"%"PF_TID_T);
        start = cnedges;

        vtx_type * pcadjncy = mycadjncy - mycxadj[mycchunkofst[*pchunkcnt-1]];
        wgt_type * pcadjwgt = mycadjwgt - mycxadj[mycchunkofst[*pchunkcnt-1]];

        int ttt = 0;
        do {
          ++ttt;
          DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],cdist),"%"PF_VTX_T);
          vtx_type const * const padjncy = (ttt==1 ? local_adjncy - cc1adjstart : gadjncy[o] - gxadj[o][gchunkofst[o][cc2]]);
          wgt_type const * const padjwgt = (ttt==1 ? local_adjwgt - cc1adjstart : gadjwgt[o] - gxadj[o][gchunkofst[o][cc2]]);

          /* transfer over vertex stuff from v and u */
          mycvwgt[c] += graph->uniformvwgt ? 1 : gvwgt[o][v];

          // (o,v) -j-> (t,k)
          // explore other coarse vertices that current coarse vertex is connected to
          for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
            k = padjncy[j];
            if (k < graph->mynvtxs[o]) {
              t = o;
            } else {
              t = gvtx_to_tid(k,graph->dist);
              k = gvtx_to_lvtx(k,graph->dist);
            }
            k = gcmap[t][k];      // now k is coarse vertex that j leads to
            if (gvtx_to_tid(k,cdist) == myid) {
              k = gvtx_to_lvtx(k,cdist);
            }
            if (k == c || k == cg) {
              /* internal edge */
            } else {
              /* external edge */
              l = k&MASK;
              i = htable[l];
              ewgt = graph->uniformadjwgt ? 1 : padjwgt[j];
              if (i == NULL_OFFSET) {
                /* new edge */
                // TODO because this is filled sequentially, we can just call
                // fwrite without manually buffering with `mycadjncy`
                pcadjncy[cnedges] = k;
                pcadjwgt[cnedges] = ewgt;
                adjwgt_sum += ewgt;
                htable[l] = (offset_type)(cnedges - start); 
                ++cnedges;  // NOTE: this is filled sequentially too
              } else {  // hash collision
                i += start;
                /* search for existing edge */
                for (jj=i;jj<cnedges;++jj) {
                  if (pcadjncy[jj] == k) {    // FIXME: we only need to store the edges of one cnode
                    pcadjwgt[jj] += ewgt;
                    adjwgt_sum += ewgt;
                    break;
                  }
                }
                if (jj == cnedges) {
                  /* we didn't find the edge, so add it */
                  // TODO call fwrite directly
                  pcadjncy[cnedges] = k;
                  pcadjwgt[cnedges] = ewgt;
                  adjwgt_sum += ewgt;
                  ++cnedges;
                }
              }
            }
          }

          v = gmatch[o][v];
          if (v >= graph->mynvtxs[o]) {
            o = gvtx_to_tid(v,graph->dist);
            v = gvtx_to_lvtx(v,graph->dist);
          }
        } while (!(myid == o && v == fcmap[c]));

        /* clear the htable */
        for (j = cnedges;j > mycxadj[c];) {
          --j;
          k = pcadjncy[j];
          l = (k&MASK);
          htable[l] = NULL_OFFSET;
        }

        mycxadj[c+1] = cnedges;

        /* check if the chunk length reaches maximum */
        size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt-1]];
        if (chunkadjs > ctrl->adjchunksize) {
#ifdef PRINT_LOOP_VAR
          fprintf(stderr, "%"PF_ADJ_T" edges written into %s\n", chunkadjs, fout1);
#endif
          mycchunkofst[*pchunkcnt] = c+1;
          ++*pchunkcnt;
          fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
          fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
        }

      }
    }
  }
  DL_ASSERT_EQUALS(c, mycnvtxs, "%"PF_VTX_T);

  size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt - 1]];
  if (chunkadjs > 0) {
#ifdef PRINT_LOOP_VAR
    fprintf(stderr, "%" PF_ADJ_T " edges written into %s\n", chunkadjs, fout1);
#endif
    mycchunkofst[*pchunkcnt] = mycnvtxs;
    fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
    fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
  } else {
    --*pchunkcnt;
  }

  fclose(cadjncy_dump);
  fclose(cadjwgt_dump);
  dl_free(local_adjncy);
  dl_free(local_adjwgt);
  dl_free(htable);

  cgraph->mynedges[myid] = cnedges;

  //graph_readjust_memory(cgraph,adjsize);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  par_chunk_graph_setup_twgts(cgraph, adjwgt_sum);
  
  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  if (myid == 0) {
    for (int i = 0; i < nthreads; ++i) {
      printf("thread %d: c#=%2zu; [%7zu|%9zu]: (", i, cgraph->chunkcnt[i], cgraph->mynvtxs[i], cgraph->mynedges[i]);
      for (int c = 0; c < cgraph->chunkcnt[i]; ++c) {
        printf("%" PF_VTX_T ",", cgraph->chunkofst[i][c + 1] - cgraph->chunkofst[i][c]);
      }
      // printf(")\n  edge count: (");
      // for (int c = 0; c < cgraph->chunkcnt[i]; ++c) {
      //   printf("%" PF_ADJ_T ",",
      //         cgraph->xadj[i][cgraph->chunkofst[i][c + 1]] - cgraph->xadj[i][cgraph->chunkofst[i][c]]);
      // }
      printf(")\n");
    }
    fflush(stdout);
  }
  dlthread_barrier(ctrl->comm);

  // DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void par_contract_chunk_graph(
    ctrl_type * const ctrl, 
    graph_type * const graph, 
    vtx_type const mycnvtxs, 
    vtx_type const * const * const gmatch, 
    vtx_type const * const fcmap)
{
  switch (ctrl->contype) {
    case MTMETIS_CONTYPE_CLS:
      S_par_contract_CLS_quadratic(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    case MTMETIS_CONTYPE_DENSE:
      exit(1);
      // S_par_contract_DENSE(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    case MTMETIS_CONTYPE_SORT:
      exit(1);
      // S_par_contract_SORT(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    default:
      dl_error("Unknown contraction type '%d'\n",ctrl->contype);
  }
}




#endif
