#!/bin/bash

if [ $# -ne 0 ];then
        echo "Usage: sh $0"
        exit 1
fi

echo [MSG] Start submitting vireo jobs...
while read line
do
	qsub vireo.pbs -v capture=`echo $line | cut -f1 -d','`,n_donor=`echo $line | cut -f2 -d','`,workDIR=`pwd`
done < "vireo3.list" #"vireo.list"

echo [MSG] All jobs submitted! 
