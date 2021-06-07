#!/bin/bash
# aedavids@ucsc.edu
# 5/23/21
#
# searches for STAR alignements and run HTSeq-count
#


set -x   # turn debug on
# set -x # turn debug off

#
# find the STAR alignments
# htseq-count requires sorted alignment for paired reads
#

annotationDir=/private/groups/kimlab/genomes.annotations
dataRoot=/private/groups/kimlab/kras.ipsc/data
cd $dataRoot

function  skipPath {
    #
    # find returns
    # two paths that end in unmapped. if not unmapped/unmapped skip
    #
    # bulk.data/day.5/ctrl.2/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped
    # bulk.data/day.5/ctrl.1/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/unmapped
    #
    path=$1

    ret=0 # false    
    # check if ends in unmapped
    echo $path | grep "unmapped$" > /dev/null 
    grepExitStatus=$? # 0 means found

    if [ $grepExitStatus -eq 0 ];
    then
        # check if end in unmapped/unmapped
        echo $path | grep "unmapped/unmapped$" > /dev/null 
        grepExitStatus2=$? # 0 means found

        if [ $grepExitStatus2 -ne 0 ];
        then
            ret=1 #true
            # bulk.data/day.5/ctrl.2/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped            
        fi
    fi

    echo $ret
}


unmappedType=unmapped
# unmappedType=m1
# unmappedType=m2
# unmappedType=m12
# unmappedType=d
# unmappedType=u

# find paths like */STAR*/unmapped/
#
# bulk.data/day.7/kras.2/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/unmapped
# bulk.data/day.7/kras.2/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/m1
STARList=`find {bulk.data,exo.data} -type d -name $unmappedType | grep STAR `
#STARList=`find {bulk.data,exo.data} -path $unmappedType | grep STAR `


# declare -i debugCount
# debugCount=0

for star in $STARList;
do
    printf "\n\n###############\n"
    
    # echo $debugCount
    
    # if [ $debugCount -ge 3 ];
    # then
    #    printf "EARLY EXIT\n"
    #    exit 0
    # fi      
    # debugCount=${debugCount}+1


    shouldSkip=`skipPath $star`
    printf "shouldSkip=$shouldSkip \n$star\n"
    if [ $shouldSkip -eq 1 ] ;
    then
        continue
    fi
    
    
    bam=${star}/Aligned.sortedByCoord.out.bam


    root=`echo $star | sed 's/STAR.*$//' `
    # bulk.data/day.7/kras.2/
    
    # htseqCount/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/m1
    suffix=`echo $star | sed 's/^.*\(STAR.*\)/\1/' `
    htSeqCountOutDir="${dataRoot}/${root}htseqCount.out/${suffix}"
    # htseqCount.out/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/m1

    if [ ! -d $htSeqCountOutDir ];
    then
        mkdir -p $htSeqCountOutDir
        htseq-count \
            --order=pos \
            $bam \
            ${annotationDir}/gencode.v37.annotation.gtf \
            > "${htSeqCountOutDir}/htseqCount.out"
    else
        printf "INFO $htSeqCountOut already exits \n"
    fi

done

    

