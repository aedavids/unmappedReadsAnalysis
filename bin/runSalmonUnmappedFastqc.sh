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
validFlags='d u m1 m2 m12 unmapped'
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

#unmappedType=unmapped
unmappedMate1List=`find {bulk.data,exo.data} -name mate1.${unmappedType}.unmapped.fastq`


# declare -i debugCount
# debugCount=0

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


    fastqOut="${sampleRoot}/fastqc.out/${salmonOutput}/unmapped/${unmappedType}"
    mkdir -p $fastqOut
                                 #mate1.d.unmapped_fastqc.zip
    expectedOutFile="${fastqOut}/mate1.${unmappedType}.unmapped_fastqc.zip"
    if [ ! -f  $expectedOutFile ];
    then
        fastqc --outdir $fastqOut --thread 6 $mate1 $mate2

    else
        printf "\n\n WARNING ${expectedOutFile} was skipped it already exists \n\n"
    fi

done

