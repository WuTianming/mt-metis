#!/bin/bash

TMPDIR=/data1/mt-metis-tmp
LOGDIR=/data1/mt-metis-final-exp/expResults-parmetis
DATE=$(date +"%Y%m%d_%H%M")

cd $TMPDIR && rm -f *
mkdir -p $LOGDIR

prefix=/data1/mt-metis-final-exp/datasets/dglpart/ogbn-

export LD_LIBRARY_PATH=/home/tmwu/local/lib

/data1/mt-metis-log/drop_page_cache
pkill nmon
data=papers100M
nmon -F $LOGDIR/${DATE}.nmon -s 2 -c 200000 -t
CMD="/bin/time -v -o $LOGDIR/${DATE}.\$MPI_LOCALRANKID.time \
        /home/tmwu/local/bin/pm_dglpart \
        ${prefix}${data} \
        1"
mpirun -np 8 bash -c "${CMD}" 2>&1 | tee $LOGDIR/${DATE}.log

pkill nmon
