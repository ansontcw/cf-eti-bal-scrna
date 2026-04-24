#!/bin/bash

if [ $# -ne 0 ];then
        echo "Usage: sh $0"
        exit 1
fi

echo [MSG] Start submitting vireo jobs...

qsub azimuth.pbs -v workDIR=`pwd`/data/SCEs

echo [MSG] All jobs submitted! 
