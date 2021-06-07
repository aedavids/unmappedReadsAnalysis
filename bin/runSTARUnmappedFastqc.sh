#!/bin/bash
# aedavids@ucsc.edu
# 5/11/21
#

scriptName=`basename $0`
if [ $# -ne 1 ]; then
    echo "ERROR: usage $scriptName fastq"
    echo "ERROR cli call: $@"
    echo "flag values are 'd' 'u' 'm1' 'm2' or 'm12' 'unmapped' "
    echo see https://salmon.readthedocs.io/en/latest/salmon.html?highlight=write#writeunmappednames
    exit 1
fi

unmappedType=$1

#
# make sure flag is valid
#
validFlags='d u m1 m2 m12'
echo $validFlags | grep $unmappedType > /dev/null
validExitStatus=$?
if [ $validExitStatus -ne 0 ]; then
    printf "ERROR '$unmappedType' is not valid, expected $validFlags\n"
    exit 1
fi

set -x   # turn debug on
# set -x # turn debug off

#
# get the unmapped salmon reads
#

dataRoot=/private/groups/kimlab/kras.ipsc/data
cd $dataRoot

unmappedMate1List=`find {bulk.data,exo.data} -name Unmapped.out.mate1`


# declare -i debugCount
# debugCount=0

for mate1 in $unmappedMate1List;
do
    printf "\n\n###############\n"
    
    echo $debugCount
    
    # if [ $debugCount -ge 3 ];
    # then
    #    printf "EARLY EXIT\n"
    #    exit 0
    # fi      
    # debugCount=${debugCount}+1


    mate2=`echo $mate1 | sed 's/mate1/mate2/g'`
    
    #
    # find the sample directory
    # easier to debug if we do not have ../../..
    #
    sampleRoot=`dirname $mate1`/../../..
    pushd $sampleRoot
    sampleRoot=`pwd`
    popd
    
    #
    # find the STAR output directory containing the unmapped read
    #
    starOutput=`dirname $mate1`/../..


    # we want the last part of the path
    # rev reverse string
    # starOutput=STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
    starOutput=`echo $starOutput | rev | cut -d / -f 5 | rev `


    fastqOut="${sampleRoot}/fastqc.out/${starOutput}/unmapped/${unmappedType}"
    mkdir -p $fastqOut
    # fastqc.out/STAR.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/unmapped/m12/Unmapped.out.mate1_fastqc.zip
    expectedOutFile="${fastqOut}/Unmapped.out.mate1_fastqc.zip"
    if [ ! -f  $expectedOutFile ];
    then
        fastqc --outdir $fastqOut --thread 6 ${mate1} ${mate2}

    else
        printf "\n\n WARNING ${expectedOutFile} was skipped it already exists \n\n"
    fi

done

