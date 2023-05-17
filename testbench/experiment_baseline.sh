#!/bin/bash

DATASET=ogbn-products bash ./launch_parmetis.sh
sleep 60

DATASET=ogbn-papers100M bash ./launch_parmetis.sh