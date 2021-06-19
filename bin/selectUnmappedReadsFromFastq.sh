#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/29/21
#
# selects reads from fastq that salmon was unable to map
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

AEWIP do not run if files already exist

source createTmpFile.sh
namesList=`createTmpFile`
add_on_exit rm $namesList
cut -d " " -f 1 ${unmappedNames} > $namesList

#
# selected the unmapped reads
#
seqtk subseq "${fastq}" "${namesList}"


