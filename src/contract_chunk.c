/**
 * @file contract_chunk.c
 * @brief Functions for performing contraction with disk-based adjlist chunks
 * @author wtm
 * @version 1
 */




#ifndef MTMETIS_CONTACT_CHUNK_C
#define MTMETIS_CONTACT_CHUNK_C




#define _GNU_SOURCE
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

struct fread_thread_data {
  FILE *fp;
  void *buffer;
  size_t buffer_size;
  size_t buffer_offset;
};

void *fread_thread(void *arg) {
  // it's soooo stupid! I want lambdas...

  struct fread_thread_data *data = (struct fread_thread_data *)arg;
  fseek(data->fp, data->buffer_offset, SEEK_SET);
  fread(data->buffer, 1, data->buffer_size, data->fp);
  return NULL;
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
    vtx_type const * const fcmap,
    int dense)
{
  int ncon = graph->ncon;

  adj_type cnedges, l, maxdeg, j, i, jj, start;
  tid_type o, t;
  vtx_type v, c, cg, k;
  wgt_type ewgt;
  graph_type * cgraph;
  graphdist_type cdist;
  offset_type * htable; // hash table
  adj_type * table;     // dense vector

  tid_type const myid = dlthread_get_id(ctrl->comm);
  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  size_t   const * const gchunkcnt = graph->chunkcnt;
  vtx_type const * const * const gchunkofst = (vtx_type const **)graph->chunkofst;
  adj_type const * const * const gxadj = (adj_type const * const *)graph->xadj;
  wgt_type const * const * const gvwgt = (wgt_type const * const *)graph->vwgt;

  vtx_type * * gadjncy = (vtx_type * *)graph->adjncy;   // the actual heap pointer may be redirected for rotating prefetch
  wgt_type * * gadjwgt = (wgt_type * *)graph->adjwgt;

  vtx_type *gadjncy_rot = vtx_alloc(ctrl->adjchunksize);
  wgt_type *gadjwgt_rot = wgt_alloc(ctrl->adjchunksize);

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

  if (maxdeg > MASK_MAX_DEG) {
    // turn on dense vector contraction
    dense = 1;
  }

  // int thresh = maxdeg / 4 + 1;

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
          : ctrl->adjchunksize;

  vtx_type * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(adjncy_chunksize);
  wgt_type * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(adjncy_chunksize);

  size_t   * cchunkcnt = cgraph->chunkcnt; cchunkcnt[myid] = 1;
  size_t   * pchunkcnt = cchunkcnt + myid;
  vtx_type * mycchunkofst = cgraph->chunkofst[myid] = vtx_alloc(maxchunkcnt * 1.5 + 5);
  mycchunkofst[0] = 0;

  if (dense)
    table = adj_init_alloc(NULL_ADJ, graph->gnvtxs);
  else
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

    /* if I don't have more chunks to process, make sure to mirror the barrier
     * calls as other threads do. Otherwise the lockstep might break and mess up
     * the rest of the program */
    if (cc1 >= gchunkcnt[myid]) {
      for (int cc2=cc1; cc2<maxchunkcnt; ++cc2) {
        dlthread_barrier(graph->comm);    // barrier every time before file read
        dlthread_barrier(graph->comm);    // after file read, too
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
        fprintf(stderr, "#%"PF_TID_T": this chunk is [%"PF_VTX_T", %"PF_VTX_T")\n", myid, cc1start, cc1end);
#endif

    vtx_type cc2start, cc2end;
    adj_type cc2adjstart, cc2adjend;
    cc2start = gchunkofst[myid][cc1], cc2end = gchunkofst[myid][cc1+1];
    cc2adjstart = gxadj[myid][cc2start], cc2adjend = gxadj[myid][cc2end];

    pthread_t read_thread[2];
    struct fread_thread_data args[2];
    args[0].buffer = gadjncy_rot;
    args[0].buffer_offset = sizeof(vtx_type) * cc2adjstart;
    args[0].buffer_size = sizeof(vtx_type) * (cc2adjend - cc2adjstart);
    args[0].fp = dadjncy_read;
    args[1].buffer = gadjwgt_rot;
    args[1].buffer_offset = sizeof(wgt_type) * cc2adjstart;
    args[1].buffer_size = sizeof(wgt_type) * (cc2adjend - cc2adjstart);
    args[1].fp = dadjwgt_read;
    pthread_create(read_thread+0, NULL, fread_thread, args+0);
    pthread_create(read_thread+1, NULL, fread_thread, args+1);

    for (int cc2=cc1; cc2<maxchunkcnt; ++cc2) {
      dlthread_barrier(graph->comm);    // barrier every time before file read
#ifdef PRINT_LOOP_VAR
      if (myid == 0) { fprintf(stderr, "chunk vector (%d,%d)\n", cc1, cc2); }
#endif
      if (cc2 < gchunkcnt[myid]) {
        // join the threads for fread
        pthread_join(read_thread[0], NULL);
        pthread_join(read_thread[1], NULL);

        void *tmp;
        tmp = gadjncy_rot; gadjncy_rot = gadjncy[myid]; gadjncy[myid] = tmp;
        tmp = gadjwgt_rot; gadjwgt_rot = gadjwgt[myid]; gadjwgt[myid] = tmp;
      }

      if (cc2+1 < gchunkcnt[myid]) {
        // launch new threads to prefetch
        cc2start = gchunkofst[myid][cc2+1], cc2end = gchunkofst[myid][cc2+2];
        cc2adjstart = gxadj[myid][cc2start], cc2adjend = gxadj[myid][cc2end];
        args[0].buffer = gadjncy_rot;
        args[0].buffer_offset = sizeof(vtx_type) * cc2adjstart;
        args[0].buffer_size = sizeof(vtx_type) * (cc2adjend - cc2adjstart);
        args[0].fp = dadjncy_read;
        args[1].buffer = gadjwgt_rot;
        args[1].buffer_offset = sizeof(wgt_type) * cc2adjstart;
        args[1].buffer_size = sizeof(wgt_type) * (cc2adjend - cc2adjstart);
        args[1].fp = dadjwgt_read;
        pthread_create(read_thread+0, NULL, fread_thread, args+0);
        pthread_create(read_thread+1, NULL, fread_thread, args+1);
      }

      dlthread_barrier(ctrl->comm);   // sync all threads to finish the read

      // coarse vertices have their info gathered in order
      // start where we left off at coarse idx = c
      for (;c<mycnvtxs;++c) {
        cg = lvtx_to_gvtx(c,myid,cdist);
        /* initialize the coarse vertex */
        for (t = 0; t < ncon; ++t) {
          mycvwgt[c * ncon + t] = 0;
        }
        // mycvwgt[c] = 0;

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

        /* check if the chunk length reaches maximum */
        size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt-1]];
        size_t deg1 = gxadj[o][v+1] - gxadj[o][v],
               deg2 = gxadj[t][k+1] - gxadj[t][k];
        if (chunkadjs + deg1 + deg2 >= ctrl->adjchunksize) {
#ifdef PRINT_LOOP_VAR
          fprintf(stderr, "%"PF_ADJ_T" edges written into %s\n", chunkadjs, fout1);
#endif
          mycchunkofst[*pchunkcnt] = c;
          ++*pchunkcnt;
          fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
          fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
        }

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
          // mycvwgt[c] += graph->uniformvwgt ? 1 : gvwgt[o][v];
          for (t = 0; t < ncon; ++t) {
            mycvwgt[c * ncon + t] += gvwgt[o][v * ncon + t];
          }

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
            vtx_type gk = k;
            k = gcmap[t][k];      // now k is coarse vertex that j leads to
            if (gvtx_to_tid(k,cdist) == myid) {
              k = gvtx_to_lvtx(k,cdist);
            }
            if (k == c || k == cg) {
              /* internal edge */
            } else {
              /* external edge */
              if (dense) {
                i = table[k];
              } else {
                l = k&MASK;
                i = htable[l];
              }
              ewgt = graph->uniformadjwgt ? 1 : padjwgt[j];
              if ((dense && i == NULL_ADJ) || (!dense && i == NULL_OFFSET)) {
                /* new edge */
                pcadjncy[cnedges] = k;
                pcadjwgt[cnedges] = ewgt;
                if (dense) {
                  table[k] = cnedges;
                } else {
                  htable[l] = (offset_type)(cnedges - start); 
                }
                ++cnedges;  // NOTE: this is filled sequentially too
              } else {      // hash collision or duplicate edge
                if (dense) {
                  pcadjwgt[i] += ewgt;
                } else {
                  // maybe hash collision? check edge list
                  i += start;
                  /* search for existing edge */
                  for (jj=i;jj<cnedges;++jj) {
                    if (pcadjncy[jj] == k) {
                      pcadjwgt[jj] += ewgt;
                      break;
                    }
                  }
                  if (jj == cnedges) {
                    /* we didn't find the edge, so add it */
                    pcadjncy[cnedges] = k;
                    pcadjwgt[cnedges] = ewgt;
                    ++cnedges;
                  }
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
        if (dense) {
          for (j = cnedges; j > mycxadj[c];) {
            --j;
            k = pcadjncy[j];
            table[k] = NULL_ADJ;
          }
        } else {
          for (j = cnedges;j > mycxadj[c];) {
            --j;
            k = pcadjncy[j];
            l = (k&MASK);
            htable[l] = NULL_OFFSET;
          }
        }

        mycxadj[c+1] = cnedges;

      } // end current coarse node c
    } // end cc2
  } // end cc1

  DL_ASSERT_EQUALS(c, mycnvtxs, "%"PF_VTX_T);

  // write left over segments of edge array
  size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt - 1]];
  if (chunkadjs > 0) {
#ifdef PRINT_LOOP_VAR
    // fprintf(stderr, "%" PF_ADJ_T " edges written into %s\n", chunkadjs, fout1);
#endif
    mycchunkofst[*pchunkcnt] = mycnvtxs;
    fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
    fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
  } else {
    --*pchunkcnt;
  }

  fclose(cadjncy_dump); fclose(cadjwgt_dump);
  fclose(dadjncy_read); fclose(dadjwgt_read);

  dl_free(local_adjncy);
  dl_free(local_adjwgt);
  dl_free(gadjncy_rot);
  dl_free(gadjwgt_rot);
  if (dense) {
    dl_free(table);
  } else {
    dl_free(htable);
  }

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

  // we may dump a certain level of graph in the format of binary graph input,
  // to be re-used or inspected later
#if 0
{
  // int dumplevel = 10;
  // int dumplevel = 4;
  int dumplevel = -1;

  if (myid == 0 && cgraph->level == dumplevel) {
    printf("writing contraction graph to file...\n");
    vtx_type voff, idx;
    wgt_type cut;
    adj_type *xadj;
    vtx_type *adjncy;
    wgt_type *adjwgt, *vwgt;
    pid_type *where = NULL, **r_where;

    tid_type const myid = dlthread_get_id(ctrl->comm);
    tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

    vtx_type const nvtxs = cgraph->nvtxs;
    int const ncon = cgraph->ncon;

    size_t const tcuts = ctrl->ninitsolutions;

    size_t myncuts = (tcuts / nthreads);

    if (myid == 0) {
    dl_start_timer(&ctrl->timers.initpart);
    }

    r_where = dlthread_get_shmem(sizeof(pid_type *), ctrl->comm);

    par_graph_gather(cgraph, &xadj, &adjncy, &vwgt, &adjwgt, &voff);

    char fname[256];
    FILE *f;
    sprintf(fname, "dcontract_lvl%d_meta.txt", dumplevel);
    f = fopen(fname, "w");
    adj_type edges = xadj[nvtxs];
    fprintf(f, "%"PF_VTX_T"\n%"PF_ADJ_T"\n%"PF_ADJ_T"\n%d\n", nvtxs, edges/2, edges, ncon);
    fclose(f);
    sprintf(fname, "dcontract_lvl%d_indptr.bin", dumplevel);
    f = fopen(fname, "w");
    fwrite(xadj, sizeof(adj_type), nvtxs+1, f);
    fclose(f);
    sprintf(fname, "dcontract_lvl%d_indices.bin", dumplevel);
    f = fopen(fname, "w");
    fwrite(adjncy, sizeof(vtx_type), edges, f);
    fclose(f);
    sprintf(fname, "dcontract_lvl%d_vwgt.bin", dumplevel);
    f = fopen(fname, "w");
    fwrite(vwgt, sizeof(wgt_type), nvtxs*ncon, f);
    fclose(f);

    abort();
  }
}
#endif

  if (myid == 0) {
    for (int i = 0; i < nthreads; ++i) {
      printf("thread %d: c#=%2zu; [%7zu|%9zu]: (", i, cgraph->chunkcnt[i], cgraph->mynvtxs[i], cgraph->mynedges[i]);
      for (int c = 0; c < cgraph->chunkcnt[i]; ++c) {
        printf("%" PF_VTX_T ",", cgraph->chunkofst[i][c + 1] - cgraph->chunkofst[i][c]);
      }
      printf(")\n");
    }
    fflush(stdout);
  }
  dlthread_barrier(ctrl->comm);

  // DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}


// two new schemes:
// 1. _random_read
// 2. _random_write

/**
 * the idea is: SSDs can handle concurrent random reads.
 * try to exploit random reads to reduce useless disk read
 * 
 * update: now use io_uring to do prefetch
*/

// everything in the mixin file is static to the current file
#include "uring_mixin.c"

static void S_par_contract_CLS_random_read(
    ctrl_type * const ctrl, 
    graph_type * const graph, 
    vtx_type const mycnvtxs, 
    vtx_type const * const * const gmatch, 
    vtx_type const * const fcmap,
    int dense)
{
  int ncon = graph->ncon;

  adj_type cnedges, l, maxdeg, j, i, jj, start;
  tid_type o, t;
  vtx_type v, c, cg, k;
  wgt_type ewgt;
  graph_type * cgraph;
  graphdist_type cdist;
  offset_type * htable; // hash table
  adj_type * table;     // dense vector

  tid_type const myid = dlthread_get_id(ctrl->comm);
  tid_type const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_type const mynvtxs = graph->mynvtxs[myid];
  size_t   const * const gchunkcnt = graph->chunkcnt;
  vtx_type const * const * const gchunkofst = (vtx_type const **)graph->chunkofst;
  adj_type const * const * const gxadj = (adj_type const * const *)graph->xadj;
  wgt_type const * const * const gvwgt = (wgt_type const * const *)graph->vwgt;
  vtx_type * const * const gadjncy = (vtx_type * const *)graph->adjncy;
  wgt_type * const * const gadjwgt = (wgt_type * const *)graph->adjwgt;

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

  if (maxdeg > MASK_MAX_DEG) {
    // turn on dense vector contraction
    dense = 1;
  }

  // int thresh = maxdeg / 4 + 1;

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
          : ctrl->adjchunksize;

  /* where to write the coarse graph */
  vtx_type * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(adjncy_chunksize);
  wgt_type * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(adjncy_chunksize);

  size_t   * cchunkcnt = cgraph->chunkcnt; cchunkcnt[myid] = 1;
  size_t   * pchunkcnt = cchunkcnt + myid;
  vtx_type * mycchunkofst = cgraph->chunkofst[myid] = vtx_alloc(maxchunkcnt * 1.5 + 5);
  mycchunkofst[0] = 0;

  if (dense)
    table = adj_init_alloc(NULL_ADJ, graph->gnvtxs);
  else
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

  struct pagecacher *dadjncy_uring = malloc(sizeof(struct pagecacher) * nthreads);
  struct pagecacher *dadjwgt_uring = malloc(sizeof(struct pagecacher) * nthreads);

  FILE **dadjncy_read, **dadjwgt_read;    // each thread reads from all files
  dadjncy_read = malloc(sizeof(FILE *) * nthreads);
  dadjwgt_read = malloc(sizeof(FILE *) * nthreads);
  // adj_type pagesize = 128*1024/8;
  adj_type pagesize = 16*1024/8;
  {
    char fin1[1024], fin2[1024];
    for (int id = 0; id < nthreads; ++id) {
      // TODO: need to open the file in binary & direct IO mode
      // TODO: 简单写！每个进程都打开所有文件，然后每个进程读取自己的文件
      sprintf(fin1, "dump_graph%zu_dadjncy_tid%"PF_TID_T".bin", graph->level, id);
      sprintf(fin2, "dump_graph%zu_dadjwgt_tid%"PF_TID_T".bin", graph->level, id);
      int fdadjncy = open(fin1, O_RDONLY | O_DIRECT);
      int fdadjwgt = open(fin2, O_RDONLY | O_DIRECT);
      DL_ASSERT(fdadjncy != -1, "open dadjncy file");
      DL_ASSERT(fdadjwgt != -1, "open dadjwgt file");
      adj_type capacity = maxdeg * 10 / pagesize + 1;
      if (ctrl->adjchunksize / nthreads / pagesize + 1 > capacity)
        capacity = ctrl->adjchunksize / nthreads / pagesize + 1;

      // call the contructor
      PageCacher(dadjncy_uring + id, capacity, graph->mynedges[id], pagesize, fdadjncy);
      PageCacher(dadjwgt_uring + id, capacity, graph->mynedges[id], pagesize, fdadjwgt);

      FILE *fadjncy = fopen(fin1, "rb");  // a separate file descriptor
      FILE *fadjwgt = fopen(fin2, "rb");
      DL_ASSERT(fadjncy != NULL, "open dadjncy file");
      DL_ASSERT(fadjwgt != NULL, "open dadjwgt file");
      dadjncy_read[id] = fadjncy;
      dadjwgt_read[id] = fadjwgt;
    }
  }
  dlthread_barrier(ctrl->comm);

  c = 0;
  twgt_type adjwgt_sum = 0;

  for (int cc1=0; cc1<maxchunkcnt; ++cc1) {
    if (cc1 >= gchunkcnt[myid]) {
      continue;
    }

    // read the chunk into thread-shared mem
    vtx_type cc1start = gchunkofst[myid][cc1], cc1end = gchunkofst[myid][cc1+1];
    adj_type cc1adjstart = gxadj[myid][cc1start], cc1adjend = gxadj[myid][cc1end];
    fseek(dadjncy_read[myid], sizeof(vtx_type) * cc1adjstart, SEEK_SET);
    fseek(dadjwgt_read[myid], sizeof(wgt_type) * cc1adjstart, SEEK_SET);
    fread(gadjncy[myid], sizeof(vtx_type), cc1adjend - cc1adjstart, dadjncy_read[myid]);
    fread(gadjwgt[myid], sizeof(wgt_type), cc1adjend - cc1adjstart, dadjwgt_read[myid]);

    // start where we left off at coarse idx = c
    for (;c<mycnvtxs;++c) {
      v = fcmap[c];
      if (v >= cc1end) {  // I'm done with the cc1 chunk
        break;
      }

      cg = lvtx_to_gvtx(c,myid,cdist);
      /* initialize the coarse vertex */
      for (t = 0; t < ncon; ++t) {
        mycvwgt[c * ncon + t] = 0;
      }

      o = myid;
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],cdist),"%"PF_VTX_T);

      k = gmatch[o][v]; t = myid;
      if (k >= graph->mynvtxs[o]) {
        t = gvtx_to_tid(k, graph->dist);
        k = gvtx_to_lvtx(k, graph->dist);
      } // now (o,v) is paired to (t,k)

      size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt-1]];
      size_t deg1 = gxadj[o][v+1] - gxadj[o][v],
             deg2 = gxadj[t][k+1] - gxadj[t][k];

      vtx_type oo = t, ov = k;

      DL_ASSERT_EQUALS(gcmap[o][v],gcmap[t][k],"%"PF_VTX_T);
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[t][k],cdist),"%"PF_VTX_T);
      DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist), graph->dist),"%"PF_TID_T);

      // do prefetch by pre-calculating the paired vertex for coarse nodes ahead
      int prefetch_depth = 64;
      for (int i=0; i<prefetch_depth; ++i) {
        vtx_type cc = c + i;
        if (cc >= mycnvtxs) {
          break;
        }
        vtx_type vv = fcmap[cc];
        vtx_type kk = gmatch[myid][vv];
        tid_type tt = myid;
        if (kk >= graph->mynvtxs[myid]) {
          tt = gvtx_to_tid(kk, graph->dist);
          kk = gvtx_to_lvtx(kk, graph->dist);
        } // now (myid,vv) is paired to (tt,kk)
        adj_type start = gxadj[tt][kk];
        adj_type end = gxadj[tt][kk+1];
        for (int64_t p=start/pagesize; p<=(end-1)/pagesize; p++) {
          prefetch_page(dadjncy_uring + tt, p, false);
          prefetch_page(dadjwgt_uring + tt, p, false);
        }
      }

      /* check if the output chunk length reaches maximum */
      if (chunkadjs + deg1 + deg2 >= ctrl->adjchunksize) {
        mycchunkofst[*pchunkcnt] = c;
        ++*pchunkcnt;
        fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
        fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
      }

      start = cnedges;

      vtx_type * pcadjncy = mycadjncy - mycxadj[mycchunkofst[*pchunkcnt-1]];
      wgt_type * pcadjwgt = mycadjwgt - mycxadj[mycchunkofst[*pchunkcnt-1]];

      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],cdist),"%"PF_VTX_T);
      vtx_type const * const padjncy = gadjncy[o] - cc1adjstart;
      wgt_type const * const padjwgt = gadjwgt[o] - cc1adjstart;

      /* transfer over vertex stuff from v */
      for (t = 0; t < ncon; ++t) {
        mycvwgt[c * ncon + t] += gvwgt[o][v * ncon + t];
      }

      // explore other coarse vertices that current coarse vertex is connected to
      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = padjncy[j];
        if (k < graph->mynvtxs[o]) {
          t = o;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        } // get (t,k) pair
        k = gcmap[t][k]; // use k to refer to the adjacent coarse vertex
        if (gvtx_to_tid(k,cdist) == myid) {
          k = gvtx_to_lvtx(k,cdist);
        }
        if (k == c || k == cg) {
          /* internal edge (the edge will be discarded) */
        } else {
          /* external edge */
          if (dense) {
            i = table[k];
          } else {
            l = k&MASK;
            i = htable[l];
          }
          ewgt = graph->uniformadjwgt ? 1 : padjwgt[j];
          if ((dense && i == NULL_ADJ) || (!dense && i == NULL_OFFSET)) {
            /* new edge */
            pcadjncy[cnedges] = k;
            pcadjwgt[cnedges] = ewgt;
            if (dense) {
              table[k] = cnedges;
            } else {
              htable[l] = (offset_type)(cnedges - start); 
            }
            ++cnedges;
          } else {      // hash collision or duplicate edge
            if (dense) {
              pcadjwgt[i] += ewgt;
            } else { // maybe hash collision? check edge list
              i += start;
              /* search for existing edge */
              for (jj=i;jj<cnedges;++jj) {
                if (pcadjncy[jj] == k) {
                  pcadjwgt[jj] += ewgt;
                  break;
                }
              }
              if (jj == cnedges) {
                /* we didn't find the edge, so add it */
                pcadjncy[cnedges] = k;
                pcadjwgt[cnedges] = ewgt;
                ++cnedges;
              }
            }
          }
        }
      }

      if (o == oo && v == ov) {
        // self-loop
        // do nothing
        // OPTIMIZE: pages do not have to be the same size.........
        // 可以弄一个超大page，然后后面的指针指向它的偏移就可以了
        adj_type adjstart, adjend;
        adjstart = gxadj[oo][ov], adjend = gxadj[oo][ov+1];

        for (int64_t p=adjstart/pagesize; p<=(adjend-1)/pagesize; p++) {
          discard_page(dadjncy_uring + oo, p);
          discard_page(dadjwgt_uring + oo, p);
        }
      } else {
        adj_type adjstart, adjend;
        adjstart = gxadj[oo][ov], adjend = gxadj[oo][ov+1];

        for (int64_t p=adjstart/pagesize; p<=(adjend-1)/pagesize; p++) {
          fetch_page(dadjncy_uring + oo, p);
          fetch_page(dadjwgt_uring + oo, p);
        }

        o = oo, v = ov;

        /* transfer over vertex stuff from v */
        for (t = 0; t < ncon; ++t) {
          mycvwgt[c * ncon + t] += gvwgt[o][v * ncon + t];
        }

        // below is copied, only altering the way to access the edge list
        for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
          k = dadjncy_uring[o].cached_pages[j/pagesize][j%pagesize];
          if (k < graph->mynvtxs[o]) {
            t = o;
          } else {
            t = gvtx_to_tid(k,graph->dist);
            k = gvtx_to_lvtx(k,graph->dist);
          } // get (t,k) pair
          k = gcmap[t][k]; // use k to refer to the adjacent coarse vertex
          if (gvtx_to_tid(k,cdist) == myid) {
            k = gvtx_to_lvtx(k,cdist);
          }
          if (k == c || k == cg) {
            /* internal edge (the edge will be discarded) */
          } else {
            /* external edge */
            if (dense) {
              i = table[k];
            } else {
              l = k&MASK;
              i = htable[l];
            }
            ewgt = graph->uniformadjwgt ? 1 : dadjwgt_uring[o].cached_pages[j/pagesize][j%pagesize];
            if ((dense && i == NULL_ADJ) || (!dense && i == NULL_OFFSET)) {
              /* new edge */
              pcadjncy[cnedges] = k;
              pcadjwgt[cnedges] = ewgt;
              if (dense) {
                table[k] = cnedges;
              } else {
                htable[l] = (offset_type)(cnedges - start); 
              }
              ++cnedges;
            } else {      // hash collision or duplicate edge
              if (dense) {
                pcadjwgt[i] += ewgt;
              } else { // maybe hash collision? check edge list
                i += start;
                /* search for existing edge */
                for (jj=i;jj<cnedges;++jj) {
                  if (pcadjncy[jj] == k) {
                    pcadjwgt[jj] += ewgt;
                    break;
                  }
                }
                if (jj == cnedges) {
                  /* we didn't find the edge, so add it */
                  pcadjncy[cnedges] = k;
                  pcadjwgt[cnedges] = ewgt;
                  ++cnedges;
                }
              }
            }
          }
        }

        for (int64_t p=adjstart/pagesize; p<=(adjend-1)/pagesize; p++) {
          discard_page(dadjncy_uring + oo, p);
          discard_page(dadjwgt_uring + oo, p);
        }
      }

      /* clear the htable */
      if (dense) {
        for (j = cnedges; j > mycxadj[c];) {
          --j;
          k = pcadjncy[j];
          table[k] = NULL_ADJ;
        }
      } else {
        for (j = cnedges;j > mycxadj[c];) {
          --j;
          k = pcadjncy[j];
          l = (k&MASK);
          htable[l] = NULL_OFFSET;
        }
      }

      mycxadj[c+1] = cnedges;
    } // end current coarse node c
  } // end cc1

  DL_ASSERT_EQUALS(c, mycnvtxs, "%"PF_VTX_T);

  // write left over segments of edge array
  size_t chunkadjs = cnedges - mycxadj[mycchunkofst[*pchunkcnt - 1]];
  if (chunkadjs > 0) {
#ifdef PRINT_LOOP_VAR
    // fprintf(stderr, "%" PF_ADJ_T " edges written into %s\n", chunkadjs, fout1);
#endif
    mycchunkofst[*pchunkcnt] = mycnvtxs;
    fwrite(mycadjncy, sizeof(adj_type), chunkadjs, cadjncy_dump);
    fwrite(mycadjwgt, sizeof(wgt_type), chunkadjs, cadjwgt_dump);
  } else {
    --*pchunkcnt;
  }

  fclose(cadjncy_dump); fclose(cadjwgt_dump);
  for (int i = 0; i < nthreads; i++) {
    fclose(dadjncy_read[i]); fclose(dadjwgt_read[i]);
    PageCacher_Del(dadjncy_uring + i); PageCacher_Del(dadjwgt_uring + i);
  }
  free(dadjncy_read); free(dadjwgt_read);
  free(dadjncy_uring); free(dadjwgt_uring);

  if (dense) {
    dl_free(table);
  } else {
    dl_free(htable);
  }

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
    case MTMETIS_CONTYPE_CLS_QUADRATIC:
      S_par_contract_CLS_quadratic(ctrl,graph,mycnvtxs,gmatch,fcmap, /*dense*/0);
      break;
    case MTMETIS_CONTYPE_CLS_RDREAD:
      S_par_contract_CLS_random_read(ctrl,graph,mycnvtxs,gmatch,fcmap, /*dense*/0);
      break;
    case MTMETIS_CONTYPE_CLS_RDWRITE:
      // TODO
      break;
    case MTMETIS_CONTYPE_DENSE:
      S_par_contract_CLS_quadratic(ctrl,graph,mycnvtxs,gmatch,fcmap, /*dense*/1);
      // dl_error("Unsupported contraction type '%d' for chunk-based contraction\n",ctrl->contype);
      break;
    case MTMETIS_CONTYPE_SORT:
      dl_error("Unsupported contraction type '%d' for chunk-based contraction\n",ctrl->contype);
      break;
    default:
      dl_error("Unknown contraction type '%d'\n",ctrl->contype);
  }
}




#endif
