#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/29/21
#
# selects reads from fastq that salmon was able to map
#

#set -x # turn debug on
# set + x # turn debug off

scriptName=`basename $0`
if [ $# -ne 2 ]; then
    echo "ERROR: usage $scriptName fastq aux_info/unmapped_names.txt"
    echo "ERROR cli call: $@"
    exit 1
fi

fastq=$1
unmappedNames=$2

#
# get the names of the unmapped reads
#
source createTmpFile.sh
unMappedNamesList=`createTmpFile`
add_on_exit rm $unMappedNamesList
cut -d " " -f 1 ${unmappedNames} | sort > $unMappedNamesList

#
# get the names of all the reads in the fastq file
#
allNames=`createTmpFile`
add_on_exit rm $allNames

file "${fastq}" | grep gzip > /dev/null 
isGZIP=$?

if [ $isGZIP -eq 0 ]; then
    zcat $fastq | grep @ | cut -d @ -f 2 | sort > $allNames
else
    grep @ $fastq | cut -d @ -f 2 | sort > $allNames
fi

#
# get the names of the mapped reads
#
mappedNames=`createTmpFile`
add_on_exit rm $mappedNames
diff $allNames $unMappedNamesList | grep "^<" | cut -d " " -f 2 > $mappedNames

#
# select the mapped reads
#
seqtk subseq "${fastq}" "${mappedNames}"
