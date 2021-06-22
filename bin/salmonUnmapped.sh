#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/29/21
#
# defines function runSalmon()
#
# ref:
# https://salmon.readthedocs.io/en/latest/salmon.html#writemappings
# --writeMappings=AEDWIP.sam does not do what we want
# it write the salmon index reference to a file AEDWIP.sam in SAM format
#
# https://salmon.readthedocs.io/en/latest/salmon.html#writeunmappednames
# --writeUnmappedNames 
# tell Salmon to write out the names of reads (or mates in paired-end reads) that
# do not map to the transcriptome
#

set -x # turn debug on
# set + x # turn debug off

function runSalmon() {
    # runs salmon on one sample and outputs to that directory
    salmonIndexDir="$1"
    rightReads="$2"
    leftReads="$3"
    outputDir="$4"
    
    #set -x # turn debug on
    # set +x # turn debug off

    if [[ ! -f "$outputDir"/quant.sf ]]; then

        mkdir -p "$outputDir"


        # printf "##############\n"
        # printf "warning --minAssignedFrags is set to $minNumFrags to enable test data set\n"
        #         minNumFrags=1        
        #            --minAssignedFrags=$minNumFrags \        
        # printf "##############\n"
        
        #if [[ -f "$inputDir"/output_single_end.fq.gz ]]; then

            numThr=12
            salmon quant \
                   -i $salmonIndexDir \
                   --libType A \
                   -1 "${rightReads}" \
                   -2 "${leftReads}" \
                   -p $numThr \
                   --recoverOrphans \
                   --validateMappings \
                   --gcBias \
                   --seqBias \
                   --rangeFactorizationBins 4 \
                   --writeUnmappedNames \
                   --output ${outputDir} 

            salmonRet=$?
            if [ $salmonRet -ne 0 ]; then
                echo ERROR salmon "$rightReads" returned exit status "$exitStatus"
                continue
            fi

        #fi
    else
        echo "[INFO] skipping ${outputDir}/quant.sf it already exists"
    fi
}

########## main

scriptName=`basename $0`
if [ $# -ne 4 ]; then
    echo "ERROR: usage $scriptName salmonIndexDir rightReads leftReads outputDir"
    echo "ERROR cli call: $@"
    echo "see script header for more details"
    exit 1
fi

salmonIndexDir="$1"
rightReads="$2"
leftReads="$3"
outputDir="$4"

runSalmon $salmonIndexDir $rightReads $leftReads $outputDir

# set -x # turn debug on
# # set + x # turn debug off

# kl=/private/groups/kimlab
# scratch=/scratch/aedavids
# dataDir=../data

# # these are small files for debuggin the only have 4 or 5 reads
# #f1="${dataDir}/test.kras.ipsc.bulk.day.5.ctrl.1.output_forward_paired.fq"
# #f2="${dataDir}/test.kras.ipsc.bulk.day.5.ctrl.1.output_reverse_paired.fq"
# # outputDir="./test.salmon.out"

# #
# # to test salmon how salamon reports unmapped reads
# # choose a control sample with the lowest reported mapping rate
# # kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/ctrl.1	 37.5851
# #
# replicate=kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/ctrl.1
# f1="${kl}/${replicate}/output_forward_paired.fq.gz"
# f2="${kl}/${replicate}/output_reverse_paired.fq.gz"
# outputDir="../data/${replicate}.unmapped"
# mkdir -p $outputDir

# #salmonIndexDir="${kl}/indexes/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx"
# salmonIndexDir="${scratch}/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx"


# runSalmon $salmonIndexDir $f1 $f2 $outputDir

