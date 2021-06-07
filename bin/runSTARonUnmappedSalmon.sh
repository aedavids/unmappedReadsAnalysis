#!/bin/bash
# aedavids@ucsc.edu
# 5/11/21
#

set -x   # turn debug on
# set -x # turn debug off

#
# get the unmapped salmon reads
#

dataRoot=/private/groups/kimlab/kras.ipsc/data
cd $dataRoot

#
# salmon unmapped reads where split into different files
# see salmon doc for details
# unmappedType = unmapped m1, m2, m12, d, u
#

unmappedType=unmapped
unmappedType=m1
unmappedType=m2
unmappedType=m12
unmappedType=d
unmappedType=u
unmappedMate1List=`find {bulk.data,exo.data} -name mate1.${unmappedType}.unmapped.fastq`


#declare -i debugCount
#debugCount=0

for mate1 in $unmappedMate1List;
do
    printf "\n\n###############\n"
    
    # echo $debugCount
    
    # if [ $debugCount -ge 3 ];
    # then
    #    printf "EARLY EXIT\n"
    #    exit 0
    # fi      
    # debugCount=${debugCount}+1


    mate2=`echo $mate1 | sed 's/mate1/mate2/g'`
    
    #
    # find the sample directory
    # easier to debug if we do not have ../..
    #
    sampleRoot=`dirname $mate1`/../..
    pushd $sampleRoot
    sampleRoot=`pwd`
    popd
    
    #
    # find the salmon output directory containing the unmapped read
    #
    salmonOutput=`dirname $mate1`/..
    pushd $salmonOutput
    salmonOutput=`pwd`
    popd

    # we want the last part of the path
    # rev reverse string
    salmonOutput=`echo $salmonOutput | rev | cut -d / -f 1 | rev `

    #prefix=/path/to/output/dir/prefix.

    # --outFileNamePrefix is misleading, it can include path to output directory
    # we want to specify a output directory not a prefix
    starOut="${sampleRoot}/STAR.${salmonOutput}/unmapped/${unmappedType}/"
    index=/private/groups/kimlab/indexes/gen.38.p.13.v37.annotation.STAR.2.7.9a.index

    if [ ! -f "${starOut}Aligned.sortedByCoord.out.bam" ];
    then
        # create sorted by coordinate bam files required for many
        # downstream applications
        #
        # --outReadsUnmapped Fastx
        #   output in separate fasta/fastq files, Unmapped.out.mate1/2
        STAR --runThreadN 10 \
             --genomeDir $index \
             --outFileNamePrefix $starOut  \
             --outReadsUnmapped Fastx\
             --readFilesIn $mate1 $mate2 \
             --outSAMtype BAM SortedByCoordinate

        samtools index \
                 --threads 9 \
                 "${starOut}Aligned.sortedByCoord.out.bam" \
                 "${starOut}Aligned.sortedByCoord.out.bam.bai"
    else
        printf "\n\n WARNING ${starOut}Aligned.out.sam was skipped it already exists \n\n"
    fi

done

