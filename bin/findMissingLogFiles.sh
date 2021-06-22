#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
#
# 6/22/21

scriptName=`basename $0`
if [ ! $# -eq 1 ]; then
    printf "ERROR: usage $scriptName listOfSalmonLogs \n\n"
    exit 1
fi

listOfSalmonLogs=$1
if [[ ! -f $listOfSalmonLogs ]]
then
    printf "ERROR: listOfSalmonLogs $1 does not exist\n"
    exit 1
fi


salmonLogs=`cat $1`
# set -x # turn debug on
# set + x # turn debug off

exitStatus=0
for log in $salmonLogs;
do
    if [[ ! -f $log ]]
    then
        printf "[ERROR] missing $log\n"
        exitStatus=1
    fi
done

exit $exitStatus
