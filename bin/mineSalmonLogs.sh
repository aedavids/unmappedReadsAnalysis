#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/29/21
#

scriptName=`basename $0`
if [ $# -ne 1 ]; then
    echo "ERROR: usage $scriptName listOfSalmonLogs"
    echo "example of how how to create list of SalmonOutDir"
    echo " find /private/groups/kimlab/panc.plasma.2020 -name salmon_quant.log \ "
    echo "  | grep sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx \ "
    echo "  > panc.plasma.2020.sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.salmon.logs.txt "
    echo ""
    echo "tsv will be written to stdout "
    exit 1
fi

salmonLogs=`cat $1`
# set -x # turn debug on
# set + x # turn debug off

source createTmpFile.sh

# kl=/private/groups/kimlab
# krasDir="${kl}/kras.ipsc/data"
# replicates="bulk.data exo.data"
# salmonLogs=`find ${krasDir}/{"bulk.data","exo.data"} -name salmon_quant.log`

# salmonLogs="/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.te.locus.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/locus.te.combined.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/gencode.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/gencode.te.locus.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/locus.te.combined.salmon.out/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log"

#salmonLogs="/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log"

# salmonLogs=`echo $salmonLogs | cut -d " " -f 1-10`
# printf "\n\n############################\n"
# echo $salmonLogs
# printf "############################\n"


# TODO AEDWIP pick mapping rate our of aux_info/meta_info.json instead of log file
grepOutTmp=`createTmpFile`
add_on_exit rm  $grepOutTmp
grep -i "mapping rate" ${salmonLogs} > $grepOutTmp


sampleNameTmp=`createTmpFile`
printf "sampleName\n" > $sampleNameTmp
add_on_exit rm $sampleNameTmp
cut -d / $grepOutTmp -f 5,6,7,8,9 >> $sampleNameTmp

salmonOutTmp=`createTmpFile`
add_on_exit rm $salmonOutTmp
printf "salmonOut\n" > $salmonOutTmp
cut -d / $grepOutTmp -f 10 >> $salmonOutTmp

# https://stackoverflow.com/a/41996668/4586180
#echo $foo | sed -n -e 's/^.*\(\(Mapping\).*\)/\1/p'
#Mapping rate = 36.8215%

# cut is a hack to pick out the numeric value
# use tr to get rid of percent
mappingRateTmp=`createTmpFile`
add_on_exit rm $mappingRateTmp
printf "mappingRate\n" > $mappingRateTmp
sed -n -e 's/^.*\(\(Mapping\).*\)/\1/p' $grepOutTmp  | cut -d = -f2 | tr "%" " " >> $mappingRateTmp

#
# pick out the ref index and fastq files used
# we use a for loop so taht we find the correct cmd_info.json file
#
indexTmp=`createTmpFile`
add_on_exit rm $indexTmp
printf "index\n" >> $indexTmp

mates1Tmp=`createTmpFile`
add_on_exit rm $mates1Tmp
printf "mate1\n" >> $mates1Tmp

mates2Tmp=`createTmpFile`
add_on_exit rm $mates2Tmp
printf "mate2\n" >> $mates2Tmp

for i in $salmonLogs;
do
    d=`dirname ${i}`
    cmdInfoJson="${d}/../cmd_info.json"

    #
    # get the ref index
    #
    indexStr=`grep index ${cmdInfoJson}`
    # "index": "/public/groups/kimlab/indexes/gencode.32.v.1.index/",
    index=`echo $indexStr | cut -d : -f 2 | tr ",\"" " "`
    # strip of last '/' if it exists
    index=`echo $index | sed 's/\/$//'`
    echo $index >> $indexTmp

    #
    # get the fastq files
    #
    mates1Str=`grep mates1 ${cmdInfoJson}`
    mates2Str=`grep mates2 ${cmdInfoJson}`
    #     "mates2": "output_reverse_paired.fq.gz",
    mates1=`echo $mates1Str | cut -d : -f 2 | tr ",\"" " "`
    mates2=`echo $mates2Str | cut -d : -f 2 | tr ",\"" " "`    
    echo $mates1 >> $mates1Tmp
    echo $mates2 >> $mates2Tmp    

done

# printf "\n\n################### indexTmp\n"
# cat $indexTmp
# printf "\n\n################### indexTmp\n"


#
# count the number of unmapped reads
#
readCountTmp=`createTmpFile`
add_on_exit rm $readCountTmp
masterCountReads.sh $salmonLogs > $readCountTmp

#
# create output file
#
pasteTmp=`createTmpFile`
add_on_exit rm $pasteTmp

paste $sampleNameTmp $mappingRateTmp $salmonOutTmp $readCountTmp $indexTmp $mates1Tmp $mates2Tmp > $pasteTmp
cat $pasteTmp



