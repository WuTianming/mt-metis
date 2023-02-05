#ifndef MTMETIS_ASYNCIO_H

#include "partition.h"
#include "coarsen.h"
#include "initpart.h"
#include "uncoarsen.h"
#include "refine.h"
#include "check.h"
#include "imetis.h"

#include <pthread.h>

typedef struct asyncio_task {
    ctrl_type *ctrl;
    graph_type *graph;      // graph to be dumped into disk
} asyncio_task;



void async_dump_to_disk  (ctrl_type *ctrl, graph_type *graph);
void async_read_from_disk(ctrl_type *ctrl, graph_type *graph);
void await_read_from_disk(ctrl_type *ctrl, graph_type *graph);


void async_rmpath(char *fname);
void async_cmd(char *line);



#define MTMETIS_ASYNCIO_H
#endif