/**
 * @file size.c
 * @brief Example of using WildRiver to open and read a matrix.
 * @author Dominique LaSalle <dominique@solidlake.com>
 * Copyright 2017
 * @version 1
 * @date 2017-07-25
 */




#include <stdio.h>
#include "wildriver.h"




int main(
  int argc,
  char ** argv)
{
  int rc;
  wildriver_matrix_handle * handle;  
  char * filename;
  size_t * rowptr;
  uint32_t * ind;
  double * val; 

  rc = 0;
  handle = NULL;

  if (argc != 2) {
    printf("Must supply a single matrix/graph file to read.\n");
    rc = 1;
    goto CLEANUP;
  }

  filename = argv[1];

  if ((handle = wildriver_open_matrix(filename, WILDRIVER_IN)) == NULL) {
    printf("Failed to open '%s'.", filename);
    rc = 1;
    goto CLEANUP;
  }

  if ((rowptr = malloc(sizeof(*rowptr)*(handle->nrows+1))) == NULL) {
    printf("Failed to allocate 'rowptr'.");
    rc = 1;
    goto CLEANUP;
  }

  if ((ind = malloc(sizeof(*ind)*(handle->nnz))) == NULL) {
    printf("Failed to allocate 'ind'.");
    rc = 1;
    goto CLEANUP;
  }

  if ((val = malloc(sizeof(*val)*(handle->nnz))) == NULL) {
    printf("Failed to allocate 'val'.");
    rc = 1;
    goto CLEANUP;
  }

  if (wildriver_load_matrix(handle, rowptr, ind, val, NULL)) {
    printf("Successfully loaded matrix of size %lux%lu with %zu " \
        "non-zeroes.\n", handle->nrows, handle->ncols, handle->nnz);
  } else {
    printf("Failed to load matrix\n");
  }

  CLEANUP:

  if (handle != NULL) {
    wildriver_close_matrix(handle);
  }
  
  return rc; 
}
