#!/bin/bash

###################### CONFIGURATIONS ########################

# where the log files are stored
LOGDIR=/data1/mt-metis-final-exp/expResults-parmetis

# where the executable is
DGLPART=/path/to/prefix/bin/pm_dglpart

# where the shared libraries are
export LD_LIBRARY_PATH=/path/to/prefix/lib:$LD_LIBRARY_PATH

# the dataset. MUST be under the current working directory
# value passed using the environment variable
data=`printenv DATASET`

# An environment variable is used to output separate results for each process.
# For Intel MPI, the environment variable `MPI_LOCALRANKID` is set.
# For OpenMPI, please use `OMPI_COMM_WORLD_LOCAL_RANK` instead.
MPI_RANK_ENV=MPI_LOCALRANKID
# MPI_RANK_ENV=OMPI_COMM_WORLD_LOCAL_RANK

##############################################################




mkdir -p $LOGDIR

DATE=$(date +"%Y%m%d_%H%M")
CMD="/bin/time -v -o $LOGDIR/${DATE}.\$${MPI_RANK_ENV}.time \
        $DGLPART \
        ${data} \
        1"
echo "# COMMANDLINE: " $CMD > $LOGDIR/${DATE}.txt
mpirun -np 8 bash -c "${CMD}" 2>&1 | tee -a $LOGDIR/${DATE}.log