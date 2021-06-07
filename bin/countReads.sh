#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/13/21
#

scriptName=$0
if [ $# -ne 1 ];
then
    printf "ERROR missing fastq file argument \n"
    exit 1
fi

fastqArg=$1

if [ ! -s ${fastqArg} ];
then
    
   # file is empty
   printf "0\n"
   exit 0
fi   

declare -i count

# is file in gzip format?
file ${fastqArg} | grep gzip > /dev/null 
grepExitStatus=$? # zero means found

if [ $grepExitStatus -eq 0 ]; then
    count=`zcat ${fastqArg} | wc -l | cut -d " " -f 1`/4
else
    count=`wc -l ${fastqArg} | cut -d " " -f 1`/4
fi

printf "$count\n"


