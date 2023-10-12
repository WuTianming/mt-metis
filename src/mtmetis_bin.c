/**
 * @file mtmetis_bin.c
 * @brief Main driver function
 * @date 2023-10-10
 * 
 * This version is trimmed down to only do post-partition refinement.
 */





#ifndef MTMETIS_BIN_C
#define MTMETIS_BIN_C



#define _DEFAULT_SOURCE

#include "base.h"
#include "partition.h"
#include "order.h"
#include "graph.h"
#include "ctrl.h"

#include <wildriver.h>

#include <fcntl.h>
#include <sys/mman.h>   // for mmap()


/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define S_ARRAY_SIZE(a) \
  (sizeof(a) > 0 ? (sizeof(a) / sizeof((a)[0])) : 0)




/******************************************************************************
* OPTIONS *********************************************************************
******************************************************************************/


static const cmd_opt_pair_t RTYPE_CHOICES[] = {
  {MTMETIS_STR_RTYPE_GREEDY,"Greedy refinement",MTMETIS_RTYPE_GREEDY}
};


static const cmd_opt_pair_t PTYPE_CHOICES[] = {
  {MTMETIS_STR_PTYPE_KWAY,"K-Way Edgecut",MTMETIS_PTYPE_KWAY},
  {MTMETIS_STR_PTYPE_ESEP,"Edge Separator",MTMETIS_PTYPE_ESEP},
  {MTMETIS_STR_PTYPE_RB,"RB Edgecut",MTMETIS_PTYPE_RB},
  {MTMETIS_STR_PTYPE_VSEP,"Vertex Separator",MTMETIS_PTYPE_VSEP},
  {MTMETIS_STR_PTYPE_ND,"Nested Dissection",MTMETIS_PTYPE_ND}
};


static const cmd_opt_pair_t VERBOSITY_CHOICES[] = {
  {MTMETIS_STR_VERBOSITY_NONE,"Do not print any runtime information.", \
      MTMETIS_VERBOSITY_NONE},
  {MTMETIS_STR_VERBOSITY_LOW,"Print only metric summary.", \
      MTMETIS_VERBOSITY_LOW},
  {MTMETIS_STR_VERBOSITY_MEDIUM,"Print summary run information.", \
      MTMETIS_VERBOSITY_MEDIUM},
  {MTMETIS_STR_VERBOSITY_HIGH,"Print verbose run information.", \
      MTMETIS_VERBOSITY_HIGH},
  {MTMETIS_STR_VERBOSITY_MAXIMUM,"Print everything.", \
      MTMETIS_VERBOSITY_MAXIMUM}
};


static const cmd_opt_pair_t DISTRIBUTION_CHOICES[] = {
  {MTMETIS_STR_DISTRIBUTION_BLOCK,"Distribute the vertices in continuous " \
      "from the initial ordering.",MTMETIS_DISTRIBUTION_BLOCK},
  {MTMETIS_STR_DISTRIBUTION_CYCLIC,"Distribute the vertices in a cyclic " \
      "fashion.",MTMETIS_DISTRIBUTION_CYCLIC},
  {MTMETIS_STR_DISTRIBUTION_BLOCKCYCLIC,"Distribute the vertices in a " \
      "blockcyclic fashion.",MTMETIS_DISTRIBUTION_BLOCKCYCLIC}
};


static const cmd_opt_pair_t IGNOREWEIGHTS_CHOICES[] = {
  {MTMETIS_STR_IGNORE_NONE,"Use all weights normally", \
      MTMETIS_IGNORE_NONE},
  {MTMETIS_STR_IGNORE_VERTEXWEIGHTS,"Force all vertex weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS},
  {MTMETIS_STR_IGNORE_EDGEWEIGHTS,"Force all edge weights to be one", \
      MTMETIS_IGNORE_EDGEWEIGHTS},
  {MTMETIS_STR_IGNORE_BOTH,"Force all weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS | MTMETIS_IGNORE_EDGEWEIGHTS}
};



static const cmd_opt_t OPTS[] = {
  {MTMETIS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_RTYPE,'r',"rtype","The type of refinement " \
      "(default=greedy).",CMD_OPT_CHOICE,RTYPE_CHOICES, \
      S_ARRAY_SIZE(RTYPE_CHOICES)},
  {MTMETIS_OPTION_SEED,'s',"seed","The random seed to use.",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NCUTS,'N',"cuts","The number of cuts to " \
      "generate using successive random seeds (default=1).",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NRUNS,'n',"runs","The number of partitionings to " \
      "generate using successive random seeds (default=1).",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NINITSOLUTIONS,'i',"initialcuts","The number of " \
      "initial cuts to generate at the coarsest level (default=8).", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_NITER,'R',"nrefpass","The maximum number of refinement " \
      "passes (default=8).",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_TIME,'t',"times","Print timing information",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_NTHREADS,'T',"threads","The number of threads to use.", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_VERBOSITY,'v',"verbosity","The amount of information to " \
      "print during partitioning (default=none).",CMD_OPT_CHOICE, \
      VERBOSITY_CHOICES,S_ARRAY_SIZE(VERBOSITY_CHOICES)},
  {MTMETIS_OPTION_DISTRIBUTION,'D',"distribution","The distribution to use " \
      "for assigning vertices to threads (default=blockcyclic).", \
      CMD_OPT_CHOICE,DISTRIBUTION_CHOICES,S_ARRAY_SIZE(DISTRIBUTION_CHOICES)},
  {MTMETIS_OPTION_UBFACTOR,'b',"balance","The balance constraint " \
      "(default=1.03, which means allowing for a 3%% imbalance).", \
      CMD_OPT_FLOAT,NULL,0},
  {MTMETIS_OPTION_ONDISK,'o',"ondisk","Store graphs to disk " \
      "during coarsening in order to reduce memory requirements (default=true).", \
      CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_ADJCHUNKSIZE,'K',"chunk","Chunk size (in Mega-edges). Only one " \
      "chunk of all edges at a time resident in memory (default=unlimited).", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_PTYPE,'p',"ptype","The type of partition to compute " \
      "(default=kway)",CMD_OPT_CHOICE,PTYPE_CHOICES, \
      S_ARRAY_SIZE(PTYPE_CHOICES)},
  {MTMETIS_OPTION_RUNSTATS,'C',"partstats","Print statics on quality of " \
      "partitions.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_METIS,'M',"metis","When run with one thread, call Metis " \
      "directly.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_LEAFMATCH,'L',"leafmatch","Match leaf vertices together " \
      "if there are too many unmatched vertices (default=true).", \
      CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_REMOVEISLANDS,'I',"removeislands","Remove island vertices " \
      "before partitioning (default=false).",CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_VWGTDEGREE,'V',"vwgtdegree","Use the degree of each " \
      "vertex as an *extra* constraint (default=false).",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_VWGTONE,'1',"vwgtone","Use 1 as an *extra* vertex weight " \
      "(default=false).",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_IGNORE,'W',"ignoreweights","Ignore input weights " \
      "on a graph file (default=none).",CMD_OPT_CHOICE,IGNOREWEIGHTS_CHOICES, \
      S_ARRAY_SIZE(IGNOREWEIGHTS_CHOICES)},
  {MTMETIS_OPTION_TRAINSETLIST,'l',"trainsetlist","Uses a file as a list " \
      "of training-set nodes. These nodes will be prepended a weight of 1.", \
      CMD_OPT_STRING,NULL,0},
  {MTMETIS_OPTION_HILLSIZE,'H',"hillsize","The limit to use when searching " \
      "for hills (default=16). This only applies to hill climbing " \
      "refinement types.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_VERSION,'\0',"version","Display the current version.", \
      CMD_OPT_FLAG,NULL,0}
};


static const size_t NOPTS = S_ARRAY_SIZE(OPTS);


#undef S_ARRAY_SIZE




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int S_usage(
    char const * const name,
    FILE * const fout)
{
  fprintf(fout,"USAGE:\n");
  fprintf(fout,"%s [options] <graphfile> <nparts> [ <partfile> | - ]\n", \
      name);
  fprintf(fout,"%s -p nd [options] <graphfile> [ <permfile> | - ]\n", \
      name);
  fprintf(fout,"\n");
  fprintf(fout,"Options:\n");
  fprint_cmd_opts(fout,OPTS,NOPTS);

  return 1;
}


static double * S_parse_args(
    cmd_arg_t * args, 
    size_t nargs,
    char const ** r_input, 
    char const ** r_output)
{
  size_t i, xarg;
  double * options = NULL;
  const char * input_file = NULL, * output_file = NULL;

  options = mtmetis_init_options();

  /* set default verbosity to low */
  options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_LOW;
  options[MTMETIS_OPTION_NPARTS] = 2.0;
  options[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_KWAY;

  for (i=0;i<nargs;++i) {
    switch (args[i].type) {
      case CMD_OPT_CHOICE:
        options[args[i].id] = (double)args[i].val.o;
        break;
      case CMD_OPT_BOOL:
        options[args[i].id] = (double)args[i].val.b;
        break;
      case CMD_OPT_INT:
        options[args[i].id] = (double)args[i].val.i;
        break;
      case CMD_OPT_FLOAT:
        options[args[i].id] = (double)args[i].val.f;
        break;
      case CMD_OPT_FLAG:
        options[args[i].id] = 1.0;
        break;
      default:
        break;
    }
  }

  xarg = 0;
  for (i=0;i<nargs;++i) {
    /* check for help */
    if (args[i].id == MTMETIS_OPTION_HELP) {
      goto CLEANUP;
    } else if (args[i].id == MTMETIS_OPTION_VERSION) {
      printf("mt-Metis %d.%d.%d\n",MTMETIS_VER_MAJOR,MTMETIS_VER_MINOR, \
          MTMETIS_VER_SUBMINOR);
      printf("Copyright 2016, The Regents of the University of Minnesota\n");
      goto CLEANUP;
    }
    if (args[i].type == CMD_OPT_XARG) {
      if (xarg == 0) {
        input_file = args[i].val.s;
      } else {
        if (options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_KWAY || \
            options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_RB) {
          if (xarg == 1) {
            options[MTMETIS_OPTION_NPARTS] = (pid_type)atoll(args[i].val.s);
          } else if (xarg == 2) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        } else {
          if (xarg == 1) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        }
      }
      ++xarg;
    }
  }

  if (input_file == NULL) {
    eprintf("Must supply at least an input graph to partition\n");
    goto CLEANUP;
  }

  *r_output = output_file;
  *r_input = input_file;

  return options;

  CLEANUP:
  dl_free(options);
  *r_output = NULL;
  *r_input = NULL;

  return NULL;
}

/* calculate statistics to show if edge count is balanced */
/* NOTE: this shouldn't be included in a git commit */
void count_total_deg_of_part(int argc, char ** argv) {
  vtx_type nvtxs, i;
  adj_type * xadj = NULL; vtx_type * adjncy = NULL;
  wgt_type * vwgt = NULL, * adjwgt = NULL;

  char const * input_file = NULL, * input_parts = NULL;
  input_file = argv[1];
  printf("reading input graph file %s...\n", input_file);
  int rv = wildriver_read_graph(input_file,&nvtxs,NULL,NULL,NULL,&xadj,&adjncy,&vwgt,&adjwgt);
  printf("stats: nvtxs = %"PF_VTX_T"\n", nvtxs);

  input_parts = argv[2];
  printf("reading input parts file %s...\n", input_parts);
  FILE *fparts = fopen(input_parts, "r");

#define DEG(k) (xadj[(k)+1] - xadj[(k)])

  int part;
  long long total[16];
  memset(total, 0, sizeof(total));
  for (size_t i = 0; i < nvtxs; ++i) {
    fscanf(fparts, "%d", &part);
    total[part] += DEG(i);
  }

  long long mx = 0;
  double avg = 0;
  for (int i = 0; i < 8; ++i) {
    printf("%lld ", total[i]);
    if (total[i] > mx) mx = total[i];
    avg += total[i];
  }
  printf("\n");
  avg /= 8;
  printf("avg = %.1lf, max = %lld\n", avg, mx);
  printf("max/avg = %.3lf\n", mx / avg);

#undef DEG

  exit(0);
}

int read_binary_dataset(
    char const * const prefix,
    vtx_type * const r_nvtxs,
    int * const r_ncon,
    adj_type ** const r_xadj,
    vtx_type ** const r_adjncy,
    wgt_type ** const r_vwgt,
    wgt_type ** const r_adjwgt,
    int usemmap)
{
  adj_type nedges;
  FILE *meta, *f_xadj, *f_adjncy, *f_vwgt;
  int fd_vwgt, fd_xadj, fd_adjncy;

  char fname[256];

  sprintf(fname, "%s_meta.txt", prefix);
  meta = fopen(fname, "r");
  if (fscanf(meta, "%"PF_VTX_T"%*"PF_ADJ_T"%"PF_ADJ_T"%d", r_nvtxs, &nedges, r_ncon) == EOF)
    dl_error("cannot open metadata %s\n", fname);
  fclose(meta);

  printf("n=%lld m=%lld ncon=%d\n", (long long)*r_nvtxs, (long long)nedges, (int)*r_ncon);

  if (*r_ncon > 0) {
    sprintf(fname, "%s_vwgt.bin", prefix);
    size_t vwgtlen = (*r_nvtxs) * (*r_ncon);
    if (usemmap) {
      fd_vwgt = open(fname, O_RDONLY);
      *r_vwgt = mmap(NULL, sizeof(wgt_type) * vwgtlen, PROT_READ, MAP_PRIVATE, fd_vwgt, 0);
      close(fd_vwgt);
    } else {
      *r_vwgt = wgt_alloc(vwgtlen);
      f_vwgt = fopen(fname, "rb");
      if (fread(*r_xadj, sizeof(adj_type), vwgtlen, f_vwgt) != vwgtlen)
        dl_error("cannot read file %s for vertex weight\n", fname);
      fclose(f_vwgt);
    }
  } else {
    *r_vwgt = NULL;
  }

  sprintf(fname, "%s_indptr.bin", prefix);
  size_t xadjlen = (*r_nvtxs + 1);
  if (usemmap) {
    fd_xadj = open(fname, O_RDONLY);
    *r_xadj = mmap(NULL, sizeof(adj_type) * xadjlen, PROT_READ, MAP_PRIVATE, fd_xadj, 0);
    close(fd_xadj);
  } else {
    *r_xadj = adj_alloc(xadjlen);
    f_xadj = fopen(fname, "rb");
    if (fread(*r_xadj, sizeof(adj_type), xadjlen, f_xadj) != xadjlen)
      dl_error("cannot read file %s for CSR indptr\n", fname);
    fclose(f_xadj);
  }

  // madvise(*r_xadj, xadjlen, );

  sprintf(fname, "%s_indices.bin", prefix);
  size_t adjncylen = nedges;
  if (usemmap) {
    fd_adjncy = open(fname, O_RDONLY);
    *r_adjncy = mmap(NULL, sizeof(vtx_type) * adjncylen, PROT_READ, MAP_PRIVATE, fd_adjncy, 0);
    close(fd_adjncy);
  } else {
    *r_adjncy = vtx_alloc(nedges);
    f_adjncy = fopen(fname, "rb");
    if (fread(*r_adjncy, sizeof(vtx_type), adjncylen, f_adjncy) != adjncylen)
      dl_error("cannot read file %s for CSR indices\n", fname);
    fclose(f_adjncy);
  }

  *r_adjwgt = NULL;

  return 1;
}


/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


int main(
    int argc, 
    char ** argv) 
{
  // count_total_deg_of_part(argc, argv);

  int rv, times, verbosity;
  size_t nargs;
  vtx_type nvtxs, i;
  int ncon = 1;
  adj_type * xadj = NULL;
  vtx_type * adjncy = NULL;
  wgt_type * vwgt = NULL, * adjwgt = NULL;
  double * options = NULL;
  cmd_arg_t * args = NULL;
  pid_type * owhere = NULL;
  char const * output_file = NULL, * input_file = NULL;
  dl_timer_t timer_input, timer_output;

  /* parse user specified options */
  rv = cmd_parse_args(argc-1,argv+1,OPTS,NOPTS,&args,&nargs);
  if (rv != DL_CMDLINE_SUCCESS) {
    S_usage(argv[0],stderr);
    rv = 1;
    goto CLEANUP;
  }
  options = S_parse_args(args,nargs,&input_file,&output_file);
  if (options == NULL) {
    S_usage(argv[0],stderr);
    rv = 2;
    goto CLEANUP;
  }

  /* parse verbosity and timing */
  times = options[MTMETIS_OPTION_TIME];
  verbosity = options[MTMETIS_OPTION_VERBOSITY];

  dl_init_timer(&timer_input);
  dl_init_timer(&timer_output);

  /* start timers */
  dl_start_timer(&timer_input);

  /* read the input graph */
  int load_binary = (strstr(input_file, ".graph") == NULL); // binary, if not .graph suffixed
  // int load_binary = 1;
  int is_mmaped = load_binary ? 1 : 0;

  if (load_binary) {
    // initiate xadj, adjncy, vwgt and adjwgt
    printf("Loading binary dataset...\n");
    rv = read_binary_dataset(input_file, &nvtxs, &ncon, &xadj, &adjncy, &vwgt, &adjwgt, is_mmaped);
  } else {
    /**
     * the vwgt array (vertex weights) is laid out as follows (example ncon=3):
     * 0 0 0 1 1 1 2 2 2 3 3 3
     * where "0 0 0" is the three weights/constraints of vertex 0.
     * 
     * the wildriver_read_graph function will allocate the vwgt array for us, and set ncon according to the input file.
     */
    rv = wildriver_read_graph(input_file,&nvtxs,NULL,&ncon,NULL,&xadj,&adjncy, \
        &vwgt,&adjwgt);
  }
  
  /* NOTE about input structures
    arrays from the input are:
    - xadj
    - adjncy
    - vwgt
    - adjwgt
  */

  if (rv != 1) {
    rv = 4;
    goto CLEANUP;
  }

  if (options[MTMETIS_OPTION_IGNORE] != MTMETIS_VAL_OFF) {
    if (((int)options[MTMETIS_OPTION_IGNORE]) & MTMETIS_IGNORE_EDGEWEIGHTS) {
      if (adjwgt) {
        if (!is_mmaped) dl_free(adjwgt);
        adjwgt = NULL;
      }
    }
    if (((int)options[MTMETIS_OPTION_IGNORE]) & MTMETIS_IGNORE_VERTEXWEIGHTS) {
      if (vwgt) {
        if (!is_mmaped) dl_free(vwgt);
        else munmap(vwgt, sizeof(wgt_type) * nvtxs * ncon);
        vwgt = NULL;
      }
    }
  }

  if (options[MTMETIS_OPTION_VWGTDEGREE] != MTMETIS_VAL_OFF) {
    // the degree of each node is appended as an *extra* constraint
    wgt_type * new_vwgt = wgt_alloc(nvtxs * (ncon+1));
    for (i=0;i<nvtxs;++i) {
      for (int j=0;j<ncon;++j) {
        new_vwgt[i*(ncon+1)+j] = vwgt[i*ncon+j];
      }
      new_vwgt[i*(ncon+1)+ncon] = xadj[i+1] - xadj[i];
    }
    if (vwgt) if (!is_mmaped) free(vwgt);
    vwgt = new_vwgt;
    ncon++;
  }

  if (options[MTMETIS_OPTION_VWGTONE] != MTMETIS_VAL_OFF) {
    // 1 (node count) is appended as an *extra* constraint
    wgt_type * new_vwgt = wgt_alloc(nvtxs * (ncon+1));
    for (i=0;i<nvtxs;++i) {
      for (int j=0;j<ncon;++j) {
        new_vwgt[i*(ncon+1)+j] = vwgt[i*ncon+j];
      }
      new_vwgt[i*(ncon+1)+ncon] = 1;
    }
    if (vwgt) if (!is_mmaped) free(vwgt);
    vwgt = new_vwgt;
    ncon++;
  }

  // check if we need to prepend weight = train-set
  for (i=0;i<nargs;++i) {
    if (args[i].type == CMD_OPT_STRING) {
      const char *filename = args[i].val.s;
      printf("using train-set: %s\n", filename);

      // prepend the weight
      {
        wgt_type *new_vwgt = wgt_alloc(nvtxs * (ncon + 1));

        for (i = 0; i < nvtxs; ++i) {
          for (int j = 0; j < ncon; ++j) {
            new_vwgt[i * (ncon + 1) + j+1] = vwgt[i * ncon + j];
          }
          new_vwgt[i * (ncon + 1) + 0] = 0;
        }

        FILE *f = fopen(filename, "r");
        if (f == NULL) {
          dl_error("cannot open training-set list: %s\n", filename);
        }
        int v;
        while (fscanf(f, "%d", &v) != EOF) {
          new_vwgt[v * (ncon + 1) + 0] = 1;
        }
        fclose(f);

        if (vwgt)
          if (!is_mmaped) free(vwgt);
        vwgt = new_vwgt;
        ncon++;
      }

      break;
    }
  }

  // after all the modification on constraints, we still have no constraints
  if (ncon == 0) {
    ncon = 1;
    // remove everything so that the code below will treat it as uniform vertex
    // weights
    if (vwgt) { if (!is_mmaped) free(vwgt); vwgt = NULL; }
    // vwgt = wgt_init_alloc(1, nvtxs);
  }

  vprintf(verbosity,MTMETIS_VERBOSITY_LOW,"Read '%s' with %"PF_VTX_T \
      " vertices and %"PF_ADJ_T" edges.\n",input_file,nvtxs,xadj[nvtxs]/2);

  owhere = pid_alloc(nvtxs);

  // read owhere from accompanying file
  {
    char fname[1024];
    sprintf(fname, "%s_where.txt", input_file);
    FILE *f = fopen(fname, "r");
    if (f == NULL) {
      dl_error("cannot open owhere file: %s\n", fname);
    }
    for (i=0;i<nvtxs;++i) {
      fscanf(f,"%d",owhere+i);
    }
  }

  dl_stop_timer(&timer_input);

  ///////////////////

  if (mtmetis_do_whatever_work_explicit(nvtxs,ncon,xadj,adjncy,is_mmaped,vwgt,adjwgt,options,
      owhere,NULL) != MTMETIS_SUCCESS) {
    rv = 3;
    goto CLEANUP;
  }

  ///////////////////

  dl_start_timer(&timer_output);

  if (output_file) {
    if (strcmp(output_file,"-") == 0) {
      /* write to stdout */
      for (i=0;i<nvtxs;++i) {
        printf("%"PF_PID_T"\n",owhere[i]);
      }
    } else {
      /* save to file */
      FILE * fout = fopen(output_file,"w");
      for (i=0;i<nvtxs;++i) {
        fprintf(fout,"%"PF_PID_T"\n",owhere[i]);
      }
      fclose(fout);
    }
  }

  dl_stop_timer(&timer_output);

  if (times) {
    dl_print_header("AUXILLARY TIME",'$');
    printf("Input: %.05fs\n",dl_poll_timer(&timer_input));
    printf("Output: %.05fs\n",dl_poll_timer(&timer_output));
    dl_print_footer('$');
  }

  CLEANUP:

  if (options) {
    dl_free(options);
  }
  // Optimization: these are cleaned early after data reading completed
  /*
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (vwgt) {
    dl_free(vwgt);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }
  */
  if (owhere) {
    dl_free(owhere);
  }
  if (args) {
    dl_free(args);
  }

  return 0;
}




#endif
