#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/13/21
#

# set -x # turn debug on
# set + x # turn debug off

# source createTmpFile.sh

# kl=/private/groups/kimlab
# krasDir="${kl}/kras.ipsc/data"
# replicates="bulk.data exo.data"
# salmonLogs=`find ${krasDir}/{"bulk.data","exo.data"} -name salmon_quant.log`

# # printf "\n\n############################\n"
# # echo $salmonLogs
# # printf "############################\n"

# #
# # count the number of reads
# #

 set -x

readFileTypes='mapped unmapped d u m1 m2 m12'
echo AEDWIP create a tmp file each row holds all counts for salmonout

salmonLog=/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log
#salmonLog=/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs-missing/salmon_quant.log


countString=""

printf "\n\n#################\n"
unmappedDir=`dirname $salmonLog`/../unmapped
if [ -d $unmappedDir ];
then
    for t in $readFileTypes;
    do
        if [ $t == "mapped" ] ;
        then
            mate1="${unmappedDir}/mate1.mapped.fastq"            
        else
            mate1="${unmappedDir}/mate1.${t}.unmapped.fastq"
        fi
        
        c=`countReads.sh ${mate1}`
        #countString="${c}\t${countString}"
        if [ $countString == "" ];
        then
           countString="${c}"
        else
            countString="${countString}\t${c}"
        fi
    done
    
else
    for t in $readFileTypes;
    do        
        countString="Na\t${countString}"
    done
fi

printf "${countString}\n"

# #done


# fastqArg=/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/../unmapped/mate1.m1.unmapped.fastq

# # empty file
# fastqArg=/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/../unmapped/mate1.m12.unmapped.fastq


