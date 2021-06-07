#!/bin/bash
#
# create STAR aligner index
# aedavids@ucsc.edu
# 5/10/21
#

#
# this take a long time. run using nohub for setsid
# nohup sh ./createSTARIndex.sh 2>&1 > createSTARIndex.sh.out &

set -x # turn debug on
# set +x # turn debug off

PATH="~/unmappedReadsAnalysis/bin/:~/extraCellularRNA/bin:${PATH}"
source createTmpFile.sh
add_on_exit dataIsUpSMS.sh 6508622638@txt.att.net $0 has completed

which STAR

STAR --version

output=/private/groups/kimlab/indexes
data=/private/groups/kimlab/genomes.annotations
indexName=gen.38.p.13.v37.annotation.STAR.2.7.9a.index

if [ -f ${output}/${indexName} ];
then
    printf "ERROR ${output}/${indexName} already exists \n"
    exit 1
fi

# run the following using nohub
STAR --runMode genomeGenerate runThreadN 10 \
    --genomeDir ${output}/${indexName} \
    --genomeFastaFiles ${data}/GRCh38.p13.genome.fa \
    --sjdbGTFfile ${data}/gencode.v37.annotation.gtf 
