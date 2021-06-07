#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/6/21

# call runSalmon on all the bulk and exo fastq

set -x   # turn debug on
# set +x # turn debug off

source salmonUnmapped.sh

# kl=/private/groups/kimlab
scratch=/scratch/aedavids
refIndexName=gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
salmonIndexDir="${scratch}/sel.align.${refIndexName}"

for fwdFastq in `cat ../data/bulk.exo.foward.fq.gz.list.txt`;
do
    printf "\n"
    f1=$fwdFastq
    f2=`echo $f1 | sed -e 's/forward/reverse/g'`
    root=`dirname $fwdFastq`
    outputDir="${root}/${refIndexName}"
    salmonUnmapped.sh $salmonIndexDir $f1 $f2 $outputDir
done

    
