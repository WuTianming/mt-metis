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

  if (ctrl->ondisk <= 0)
    return;

  // only save the first 3 graphs provides acceptable memory peak
  // if (gID >= 4) { return; }

  graph->ondisk = 1;  // TODO: make this a macro

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
  graph->ondisk = 0;
}

asyncio_task *get_async_task(ctrl_type *ctrl, graph_type *graph) {
  asyncio_task *t = malloc(sizeof(asyncio_task));   // Caveat: this is never freed, for now
  t->ctrl = ctrl;
  t->graph = graph;
  t->status = STANDBY;
  return t;
}

static void launch_dump(void *p) {
  asyncio_task *t = (asyncio_task *)p;
  t->graph->io_pid = pthread_self();
  t->status = ONGOING;
  S_ser_write_to_disk(t->ctrl, t->graph);
  t->status = DONE;
}

static void launch_cmd(void *p) {
  char *L = (char *)p;
  system(L);
  free(L);
}

void async_dump_to_disk(ctrl_type *ctrl, graph_type *graph) {
  asyncio_task *t = get_async_task(ctrl, graph);

  // nesting a pthread call inside OpenMP structures is perfectly okay
  pthread_t pid;

  // create thread to run `launch_dump()`
  int ret = pthread_create(&pid, NULL, &launch_dump, (void *)t);
  if (ret) {
    fprintf(stderr, "pthread_create for async dump to disk failed.\n");
    exit(1);
  }

  // the thread will be joined later when trying to read the file
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