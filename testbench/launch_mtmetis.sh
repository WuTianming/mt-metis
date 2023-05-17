#!/bin/bash

###################### CONFIGURATIONS ########################

TMPDIR=/data1/mt-metis-tmp
LOGDIR=/data1/mt-metis-final-exp/expResults-mtmetis
MTMETIS=/path/to/mt-metis/build/Linux-x86_64/bin/mtmetis

##############################################################


ARG_NTHREAD=8

# the following arguments are passed using the environment variables
ARG_DATASET=`printenv DATASET`
ARG_CHUNKSIZE=`printenv CHUNKSIZE`

DATE=$(date +"%Y%m%d_%H%M")

FLAGS="-vhigh -t"
FLAGS="$FLAGS -K${ARG_CHUNKSIZE}"

CMD="/usr/bin/env OMP_NUM_THREADS=$ARG_NTHREAD \
    $MTMETIS $FLAGS \
    $ARG_DATASET \
    8 \
    $LOGDIR/${DATE}.partfile \
    2>&1 | tee -a $LOGDIR/${DATE}.txt"

cd $TMPDIR && rm -f *
mkdir -p $LOGDIR
echo "# COMMANDLINE: " $CMD > $LOGDIR/${DATE}.txt

/usr/bin/time -o $LOGDIR/${DATE}.time -v bash -c "${CMD}" 2>&1