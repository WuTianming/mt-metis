#!/bin/bash

./launch_parmetis1.sh

sleep 60

./launch_parmetis2.sh

sleep 60

./experiment_mtmetis.sh