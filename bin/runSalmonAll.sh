#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/6/21

# call runSalmon on all the bulk and exo fastq
#
# output:
#   If fastq files are in a directory 'foo', the salmon output will be but in sub
#   directory with the name of the base name of the refIndexPath
#
# arguments:
#   refIndexPath
#       example:/scratch/aedavids/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
#
#   fileListPath
#       example bullk.exo.foward.fq.gz.list.txt
#       $ head bulk.exo.foward.fq.gz.list.txt
#       /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/output_forward_paired.fq.gz
#       /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/output_forward_paired.fq.gz
#       /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.3/output_forward_paired.fq.gz

#   assumptions:
#       - paired reads
#       - s/foward/reverse/ i.e. replace foward with reverse to find left mate


scriptName=`basename $0`
if [ $# -ne 2 ]; then
    echo "ERROR: usage $scriptName refIndexPath fileListPath"
    echo "ERROR cli call: $@"
    echo "see script header for more details"
    exit 1
fi

refIndexPath=$1
fileListPath=$2

set -x   # turn debug on
# # set +x # turn debug off

source createTmpFile.sh

function startBatch() { 
    fileListPath=$1
    refIndexName=$2
    refIndexDir=$3

    for fwdFastq in `cat fileListPath`;
    do
        f1=$fwdFastq
        f2=`echo $f1 | sed -e 's/forward/reverse/g'`
        root=`dirname $fwdFastq`
        outputDir="${root}/${refIndexName}"
        echo AEDWIP salmonUnmapped.sh $refIndexDir $f1 $f2 $outputDir
    done
}


#
# control level of concurrency so that other users can make progress
# split fileListPath into numBatches
#
# use cat so that result is only the number of lines, not num lines and file name
#numFiles=`cat bulk.data.foward.fq.gz.list.txt | wc -l `
numFiles=`cat ${fileListPath} | wc -l `

# if numFiles / numBatch is real, you get back the quotient, ie an integer value
numBatches=6
numLinesPerSplit=$(expr $numFiles / $numBatches)

tmpDir=`createTmpDir`
add_on_exit 'rm' -rf $tmpDir

# fileListPath might be relative
d=`pwd`
cd $tmpDir
split --lines=$numLinesPerSplit "$d/$fileListPath"


#
# run batches in background
#

cd $d
extraBin=/private/home/aedavids/extraCellularRNA/bin
now=`${extraBin}/dateStamp.sh`
logDir="./${scriptName}.${now}.log"
mkdir -p $logDir

for batchFile in `ls $tmpDir`;
do
    startBatch $batchFile $refIndexName $refIndexPath 2>&1 >  "${logDir}/${batchFile}" &
done

# do
#     printf "\n"
#     f1=$fwdFastq
#     f2=`echo $f1 | sed -e 's/forward/reverse/g'`
#     root=`dirname $fwdFastq`
#     outputDir="${root}/${refIndexName}"
#     salmonUnmapped.sh $salmonIndexDir $f1 $f2 $outputDir
# done

    
#
# clean up
#
'rm' -rf $tmpDir
