#include "asyncio.h"
#include "metislib.h"
#include "string.h"

static void S_ser_write_to_disk(
    ctrl_type * const ctrl,
    graph_type * const graph) {

  static int gID = 1;

  // size_t size_before = graph_size(graph);

  vtx_type *nvtxs, ncon;
  adj_type *nedges;
  tid_type nthreads = ctrl->nthreads;
  char outfile[1024];
  FILE *fpout, *metaout;

  // only save the first 3 graphs provides acceptable memory peak
  // if (gID >= 4) { return; }
  if (graph->level >= 16) { return; }

  if (graph->gID > 0) {
    sprintf(outfile, "dump_mtmetis.%d", graph->gID);
    gk_rmpath(outfile);
  }

  graph->gID = gID++;

  sprintf(outfile, "dump_mtmetis.%d", graph->gID);
  if ((fpout = fopen(outfile, "wb")) == NULL) abort();

  nvtxs  = graph->mynvtxs;
  nedges = graph->mynedges;
  ncon   = graph->ncon;     // now supports multiple constraints!

  {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (!graph->free_adjncy) {
        dl_free(graph->adjncy[myid]);
        graph->adjncy[myid] = NULL;
      }
      if (!graph->free_adjwgt) {
        dl_free(graph->adjwgt[myid]);
        graph->adjwgt[myid] = NULL;
      }
    }
  }

  sprintf(outfile, "dump_meta.%d.txt", graph->gID);
  if ((metaout = fopen(outfile, "w")) == NULL) abort();

  // TODO: dump metadata to facilitate checkpointing
  fprintf(metaout, "%d %d %d %d %d %d %d %d %d %d\n", graph->gID, (int)nthreads,
          (int)ncon, graph->uniformadjwgt, graph->uniformvwgt,
          graph->free_adjncy, graph->free_adjwgt, graph->free_vsize,
          graph->free_vwgt, graph->free_xadj);

  for (int i = 0; i < nthreads; ++i)
    fprintf(metaout, " %"PF_VTX_T, nvtxs[i]);
  fprintf(metaout, "\n");

  for (int i = 0; i < nthreads; ++i)
    fprintf(metaout, " %"PF_ADJ_T, nedges[i]);
  fprintf(metaout, "\n");

  for (int i = 0; i < nthreads; ++i)
    fprintf(metaout, " %zu", graph->chunkcnt[i]);
  fprintf(metaout, "\n");

  for (int i = 0; i < nthreads; ++i) {
    for (int j = 0; j <= graph->chunkcnt[i]; ++j) {
      fprintf(metaout, " %zu", graph->chunkofst[i][j]);
    }
    fprintf(metaout, "\n");
  }
  fprintf(metaout, "\n");
  fclose(metaout);

  if (graph->free_xadj) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->xadj[myid], sizeof(adj_type), nvtxs[myid]+1, fpout) != (size_t)(nvtxs[myid]+1))
        abort();
      dl_free(graph->xadj[myid]);
      graph->xadj[myid] = NULL;
    }
  }
  if (graph->free_vwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->vwgt[myid], sizeof(wgt_type), nvtxs[myid]*ncon, fpout) != (size_t)(nvtxs[myid]*ncon))
        abort();
      dl_free(graph->vwgt[myid]);
      graph->vwgt[myid] = NULL;
    }
  }
  if (graph->free_adjncy) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->adjncy[myid], sizeof(vtx_type), nedges[myid], fpout) != (size_t)nedges[myid])
        abort();
      dl_free(graph->adjncy[myid]);
      graph->adjncy[myid] = NULL;
    }
  }
  if (graph->free_adjwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->adjwgt[myid], sizeof(wgt_type), nedges[myid], fpout) != (size_t)nedges[myid])
        abort();
      dl_free(graph->adjwgt[myid]);
      graph->adjwgt[myid] = NULL;
    }
  }

  fclose(fpout);

  // TODO see: S_graph_free_part @ graph.c
  ///////// BEGIN copy from graph.c
  /* free auxillery things */
  /*
  if (graph->label) {
    dl_free(graph->label[myid]);
  }
  if (graph->group) {
    dl_free(graph->group[myid]);
  }
  */
  ///////// END copy from graph.c

  // size_t size_after = graph_size(graph);

  // printf("%d -> %d\n", (int)(size_before >> 20), (int)(size_after >> 20));

  return;

ERROR:
  printf("Failed on writing %s\n", outfile);
  fclose(fpout);
  gk_rmpath(outfile);
  // graph->ondisk = 0;
}

/**
 * intended use case:
 *   graph = par_graph_create(comm);
 *   graph = recover_metadata(ctrl, graph, gid);
*/
graph_type *  S_ser_recover_graph_metadata(
    ctrl_type * const ctrl,
    graph_type * graph,
    int gID) {

  exit(1);    // this function is unfinished

  vtx_type *nvtxs, ncon;
  adj_type *nedges;
  tid_type nthreads = ctrl->nthreads;
  char infile[1024];
  FILE *metaout;

  // FIXME: refer to graph.c:3572, par_graph_distribute()
  /** NOTE:
   * need to set graph->ondisk != RESIDENT
  */

  nvtxs  = graph->mynvtxs;
  nedges = graph->mynedges;
  ncon   = graph->ncon;     // TODO: this should be read from the metadata file

  sprintf(infile, "dump_meta.%d.txt", graph->gID);
  if ((metaout = fopen(infile, "r")) == NULL) {
    dl_error("cannot open metadata file for graph recovery\n");
  }

  fscanf(metaout, "%d %"PF_TID_T" %"PF_VTX_T" %d %d %d %d %d %d %d", &graph->gID, &nthreads,
         &ncon, &graph->uniformadjwgt, &graph->uniformvwgt,
         &graph->free_adjncy, &graph->free_adjwgt, &graph->free_vsize,
         &graph->free_vwgt, &graph->free_xadj);
  
  if (nthreads != ctrl->nthreads) { abort(); }

  for (int i = 0; i < nthreads; ++i)
    fscanf(metaout, "%"PF_VTX_T, &nvtxs[i]);

  for (int i = 0; i < nthreads; ++i)
    fscanf(metaout, "%"PF_ADJ_T, &nedges[i]);

  for (int i = 0; i < nthreads; ++i)
    fscanf(metaout, "%zu", &graph->chunkcnt[i]);

  for (int i = 0; i < nthreads; ++i) {
    graph->chunkofst[i] = vtx_alloc(graph->chunkcnt[i] + 2);

    for (int j = 0; j <= graph->chunkcnt[i]; ++j) {
      fscanf(metaout, "%zu", &graph->chunkofst[i][j]);
      if (j > 0 && graph->chunkofst[i][j-1] >= graph->chunkofst[i][j]) {
        fprintf(stderr, "error! chunkofst regression.\n");
        abort();
      }
    }
  }

  graph->ondisk = OFFLOADING;

  fclose(metaout);
}

static void S_ser_read_from_disk(
    ctrl_type * const ctrl,
    graph_type * const graph) {

  if (graph->ondisk == RESIDENT) return;

  vtx_type *nvtxs, ncon;
  adj_type *nedges;
  tid_type nthreads = ctrl->nthreads;
  char infile[1024];
  FILE *fpin;

  sprintf(infile, "dump_mtmetis.%d", graph->gID);
  if ((fpin = fopen(infile, "rb")) == NULL)
    return;

  nvtxs  = graph->mynvtxs;
  nedges = graph->mynedges;
  ncon   = graph->ncon;     // now supports multiple constraints!

  if (graph->free_xadj) {
    for (int myid = 0; myid < nthreads; ++myid) {
      graph->xadj[myid] = adj_alloc(nvtxs[myid]+1);
      if (fread(graph->xadj[myid], sizeof(adj_type), nvtxs[myid]+1, fpin) != (size_t)(nvtxs[myid]+1))
        abort();
      graph->xadj[myid][0] = 0;
    }
  }
  if (graph->free_vwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      graph->vwgt[myid] = wgt_alloc(nvtxs[myid]*ncon);
      if (fread(graph->vwgt[myid], sizeof(wgt_type), nvtxs[myid]*ncon, fpin) != (size_t)(nvtxs[myid]*ncon))
        abort();
    }
  }
  if (graph->free_adjncy) {
    for (int myid = 0; myid < nthreads; ++myid) {
      graph->adjncy[myid] = vtx_alloc(nedges[myid]);
      if (fread(graph->adjncy[myid], sizeof(vtx_type), nedges[myid], fpin) != (size_t)nedges[myid])
        abort();
    }
  } else {
    for (int myid = 0; myid < nthreads; ++myid) {
      mtmetis_vtx_type *ofst = graph->chunkofst[myid];
      mtmetis_adj_type *xadj = graph->xadj[myid];
      size_t adjncy_chunksize = xadj[ofst[1]];
      for (int i = 1; i < graph->chunkcnt[myid]; ++i) {
        dl_storemax(adjncy_chunksize, xadj[ofst[i+1]] - xadj[ofst[i]]);
      }
      graph->adjncy[myid] = vtx_alloc(adjncy_chunksize);
      graph->free_adjncy = 1;   // these arrays will be freed later by `par_graph_free()`
    }
  }
  if (graph->free_adjwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      graph->adjwgt[myid] = wgt_alloc(nedges[myid]);
      if (fread(graph->adjwgt[myid], sizeof(wgt_type), nedges[myid], fpin) != (size_t)nedges[myid])
        abort();
    }
  } else {
    for (int myid = 0; myid < nthreads; ++myid) {
      mtmetis_vtx_type *ofst = graph->chunkofst[myid];
      mtmetis_adj_type *xadj = graph->xadj[myid];
      size_t adjncy_chunksize = xadj[ofst[1]];
      for (int i = 1; i < graph->chunkcnt[myid]; ++i) {
        dl_storemax(adjncy_chunksize, xadj[ofst[i+1]] - xadj[ofst[i]]);
      }
      graph->adjwgt[myid] = wgt_alloc(adjncy_chunksize);
      graph->free_adjwgt = 1;
    }
  }

  fclose(fpin);

  // delete the checkpoint files
  printf("ondisk: deleting %s\n", infile);
  // gk_rmpath(infile);
  async_rmpath(infile);   // save a few hundred milliseconds by running in background

  graph->gID    = 0;

  return;

ERROR:
  dl_error("Failed to restore graph %s from the disk\n", infile);
  fclose(fpin);
  gk_rmpath(infile);
  // graph->ondisk = 0;
}

asyncio_task *get_async_task(ctrl_type *ctrl, graph_type *graph) {
  asyncio_task *t = malloc(sizeof(asyncio_task));
  t->ctrl = ctrl;
  t->graph = graph;
  return t;
}

static void *launch_dump(void *p) {
  asyncio_task *t = (asyncio_task *)p;

  S_ser_write_to_disk(t->ctrl, t->graph);

  return p;
}

static void *launch_read(void *p) {
  asyncio_task *t = (asyncio_task *)p;

  S_ser_read_from_disk(t->ctrl, t->graph);

  return p;
}

static void launch_cmd(void *p) {
  char *L = (char *)p;
  system(L);
  free(L);
}

void async_dump_to_disk(ctrl_type *ctrl, graph_type *graph) {
  if (ctrl->ondisk == 0)
    return;

  // the state machine
  assert(graph->ondisk == RESIDENT);
  graph->ondisk = OFFLOADING;

  // nesting a pthread call inside OpenMP structures is perfectly okay
  pthread_t pid;
  asyncio_task *t = get_async_task(ctrl, graph);

  // create thread to run `launch_dump()`
  int ret = pthread_create(&pid, NULL, &launch_dump, (void *)t);
  if (ret) {
    dl_error("pthread_create for async dump to disk failed.\n");
  }

  graph->io_pid = pid;

  // the thread will be joined later when trying to read the file
}

void async_read_from_disk(ctrl_type *ctrl, graph_type *graph) {
  if (ctrl->ondisk == 0)
    return;

  asyncio_task *pre_t;

  // the state machine:
  if (graph->ondisk == RESIDENT || graph->ondisk == LOADING)
    return;

  // now that the graph had been written onto the disk,
  // wait for the write to complete
  if (graph->ondisk == OFFLOADING) {
    pthread_join(graph->io_pid, (void *)&pre_t);
    graph->ondisk = LOADING;

    free(pre_t);
  }

  pthread_t pid;
  asyncio_task *t = get_async_task(ctrl, graph);
  int ret = pthread_create(&pid, NULL, &launch_read, (void *)t);
  if (ret) {
    dl_error("pthread_create for async read from disk failed.\n");
  }

  graph->io_pid = pid;

  // the thread will be joined when trying to ensure the file is ready
}

void await_read_from_disk(ctrl_type *ctrl, graph_type *graph) {
  if (ctrl->ondisk == 0)
    return;

  asyncio_task *pre_t;

  // the state machine:
  if (graph->ondisk == LOADING) {
    pthread_join(graph->io_pid, (void *)&pre_t);
    graph->io_pid = 0;
    graph->ondisk = RESIDENT;

    free(pre_t);
  }
}

void async_rmpath(char *fname) {
  char tmp[2048];
  sprintf(tmp, "rm -r %s", fname);
  async_cmd(tmp);
}

void async_cmd(char *line) {
  int len = strlen(line);

  // create a copy of the original string to be used in a new thread
  char *new_str = (char *)malloc(sizeof(char) * (len + 1));
  strncpy(new_str, line, len + 1);

  pthread_t pid;

  int ret = pthread_create(&pid, NULL, (void *)&launch_cmd, (void *)new_str);
  if (ret) {
    dl_error("pthread_create for async command failed. requested cmd: \"%s\"\n", new_str);
  }

  // make the thread self-destruct after the job is done
  pthread_detach(pid);
}