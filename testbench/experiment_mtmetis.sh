#!/bin/bash

# CHUNKSIZE=-1  bash ./launch_experiment.sh
# CHUNKSIZE=200 bash ./launch_experiment.sh
# CHUNKSIZE=150 bash ./launch_experiment.sh

prefix=/data1/mt-metis-final-exp/datasets/mtmetis-binary/ogbn-

DATASET=${prefix}products CHUNKSIZE=100 bash ./launch_mtmetis.sh
sleep 60

DATASET=${prefix}papers100M CHUNKSIZE=100 bash ./launch_mtmetis.sh
sleep 60

DATASET=${prefix}papers400M CHUNKSIZE=100 bash ./launch_mtmetis.sh
sleep 60

DATASET=${prefix}papers400M CHUNKSIZE=400 bash ./launch_mtmetis.sh
