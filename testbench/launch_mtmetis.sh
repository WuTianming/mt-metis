#!/bin/bash

ARG_NTHREAD=8
# ARG_DATASET=400M
# ARG_DATASET=/home/tmwu/ogbn-arxiv.graph
# ARG_DATASET=/home/tmwu/ogbn-papers100M.graph
# ARG_DATASET=/home/tmwu/ogbn-products.graph
ARG_DATASET=`printenv DATASET`
ARG_CHUNKSIZE=`printenv CHUNKSIZE`

TMPDIR=/data1/mt-metis-tmp
LOGDIR=/data1/mt-metis-final-exp/expResults-mtmetis
DATE=$(date +"%Y%m%d_%H%M")

MTMETIS=/home/tmwu/projMETIS/mt-metis/build/Linux-x86_64/bin/mtmetis
FLAGS="-vhigh -t"
FLAGS="$FLAGS -K${ARG_CHUNKSIZE}"

CMD="/usr/bin/env OMP_NUM_THREADS=$ARG_NTHREAD \
      $MTMETIS $FLAGS \
        $ARG_DATASET 8 \
          $LOGDIR/${DATE}.partfile \
            2>&1 | tee -a $LOGDIR/${DATE}.txt"

cd $TMPDIR && rm -f *
mkdir -p $LOGDIR
echo "# COMMANDLINE: " $CMD > $LOGDIR/${DATE}.txt

# /data1/mt-metis-log/drop_page_cache
pkill nmon --uid $USER
nmon -F $LOGDIR/${DATE}.nmon -s 2 -c 200000 -t
/usr/bin/time -o $LOGDIR/${DATE}.time -v bash -c "${CMD}" 2>&1
pkill nmon --uid $USER
