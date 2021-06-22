#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/29/21
#
# selects reads from fastq that salmon was unable to map
#

scriptName=`basename $0`
if [ $# -ne 3 ]; then
    echo "ERROR: usage $scriptName fastq flag aux_info/unmapped_names.txt"
    echo "ERROR cli call: $@"
    echo "flag values are 'd' 'u' 'm1' 'm2' or 'm12'"
    echo see https://salmon.readthedocs.io/en/latest/salmon.html?highlight=write#writeunmappednames
    exit 1
fi

#set -x # turn debug on
# set + x # turn debug off

fastq=$1
flag=$2
unmappedNames=$3

#
# make sure flag is valid
#
validFlags='d u m1 m2 m12'
echo $validFlags | grep $flag > /dev/null
validExitStatus=$?
if [ $validExitStatus -ne 0 ]; then
    printf "ERROR '$flag' is not valid, expected $validFlags\n"
    exit 1
fi

source createTmpFile.sh
tmpNamesList=`createTmpFile`
add_on_exit rm $tmpNamesList
grep "${flag}$" ${unmappedNames} | cut -d " " -f 1  > $tmpNamesList

#
# selected the unmapped reads
#
seqtk subseq "${fastq}" "${tmpNamesList}"

'rm' $tmpNamesList


