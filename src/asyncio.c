#include "asyncio.h"
#include "string.h"

static void S_ser_write_to_disk(
    ctrl_type * const ctrl,
    graph_type * const graph) {

  static int gID = 1;

  // size_t size_before = graph_size(graph);

  vtx_type *nvtxs, *nedges, ncon;
  tid_type nthreads = ctrl->nthreads;
  char outfile[1024];
  FILE *fpout;

  // only save the first 3 graphs provides acceptable memory peak
  // if (gID >= 4) { return; }

  if (graph->gID > 0) {
    sprintf(outfile, "dump_mtmetis.%d", graph->gID);
    gk_rmpath(outfile);
  }

  graph->gID    = gID++;
  sprintf(outfile, "dump_mtmetis.%d", graph->gID);

  if ((fpout = fopen(outfile, "wb")) == NULL)
    abort();

  nvtxs  = graph->mynvtxs;
  nedges = graph->mynedges;
  ncon   = 1;  // only 1 type of constraint is supported

  {
    // FIXME: re-allocate when attempting to read back
    for (int myid = 0; myid < nthreads; ++myid) {
      dl_free(graph->adjncy[myid]);
      graph->adjncy[myid] = NULL;
      dl_free(graph->adjwgt[myid]);
      graph->adjwgt[myid] = NULL;
    }
  }

  if (graph->free_xadj) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->xadj[myid], sizeof(adj_type), nvtxs[myid]+1, fpout) != (size_t)(nvtxs[myid]+1))
        abort();
        // goto ERROR;
      dl_free(graph->xadj[myid]);
      graph->xadj[myid] = NULL;
    }
    // dl_free(graph->xadj);
    // graph->xadj = NULL;
  }
  if (graph->free_vwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->vwgt[myid], sizeof(wgt_type), nvtxs[myid]*ncon, fpout) != (size_t)(nvtxs[myid]*ncon))
        abort();
      dl_free(graph->vwgt[myid]);
      graph->vwgt[myid] = NULL;
    }
    // dl_free(graph->vwgt);
    // graph->vwgt = NULL;
  }
  if (graph->free_adjncy) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->adjncy[myid], sizeof(vtx_type), nedges[myid], fpout) != (size_t)nedges[myid])
        abort();
      dl_free(graph->adjncy[myid]);
      graph->adjncy[myid] = NULL;
    }
    // dl_free(graph->adjncy);
    // graph->adjncy = NULL;
  }
  if (graph->free_adjwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      if (fwrite(graph->adjwgt[myid], sizeof(wgt_type), nedges[myid], fpout) != (size_t)nedges[myid])
        abort();
      dl_free(graph->adjwgt[myid]);
      graph->adjwgt[myid] = NULL;
    }
    // dl_free(graph->adjwgt);
    // graph->adjwgt = NULL;
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

static void S_ser_read_from_disk(
    ctrl_type * const ctrl,
    graph_type * const graph) {

  if (graph->ondisk == RESIDENT) return;

  vtx_type *nvtxs, *nedges, ncon;
  tid_type nthreads = ctrl->nthreads;
  char infile[1024];
  FILE *fpin;

  sprintf(infile, "dump_mtmetis.%d", graph->gID);

  if ((fpin = fopen(infile, "rb")) == NULL)
    return;

  nvtxs  = graph->mynvtxs;
  nedges = graph->mynedges;
  ncon   = 1;  // only 1 type of constraint is supported

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
  }
  if (graph->free_adjwgt) {
    for (int myid = 0; myid < nthreads; ++myid) {
      graph->adjwgt[myid] = wgt_alloc(nedges[myid]);
      if (fread(graph->adjwgt[myid], sizeof(wgt_type), nedges[myid], fpin) != (size_t)nedges[myid])
        abort();
    }
  }

  fclose(fpin);
  printf("ondisk: deleting %s\n", infile);
  // gk_rmpath(infile);
  async_rmpath(infile);   // save a few hundred milliseconds

  graph->gID    = 0;

  return;

ERROR:
  printf("Failed to restore graph %s from the disk\n", infile);
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
    fprintf(stderr, "pthread_create for async dump to disk failed.\n");
    exit(1);
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
    pthread_join(graph->io_pid, &pre_t); // free(pre_t);
    graph->ondisk = LOADING;
  }

  pthread_t pid;
  asyncio_task *t = get_async_task(ctrl, graph);
  int ret = pthread_create(&pid, NULL, &launch_read, (void *)t);
  if (ret) {
    fprintf(stderr, "pthread_create for async read from disk failed.\n");
    exit(1);
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
    pthread_join(graph->io_pid, &pre_t); // free(pre_t);
    graph->io_pid = 0;
    graph->ondisk = RESIDENT;
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

  int ret = pthread_create(&pid, NULL, &launch_cmd, (void *)new_str);
  if (ret) {
    fprintf(stderr, "pthread_create for async command failed. requested cmd: \"%s\"\n", new_str);
    exit(1);
  }

  // make the thread self-destruct after the job is done
  pthread_detach(pid);
}