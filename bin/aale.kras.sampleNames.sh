#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 6/22/21
#
# create sample names for use with mineSalmonLogs.sh

scriptName=`basename $0`
usage="usage: $scriptName dataDirPath refIndexName listOfSalmonLogs  \n\
creates a list of sample names for use with  mineSalmonLogs.sh \n\
\ndataDirPath is the root of the log file path. \n\
example: /private/groups/kimlab/aale.kras/data \n\
\nref is the name of the salmon index \n\
example ref: sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx \n\
\nexample of a file path in the listOfSalmonlogs\n \
/private/groups/kimlab/aale.kras/data/bulk.rna.seq/aale/input/ctrl.1/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/logs/salmon_quant.log\n\
\n The sample will be the path string between dataDirPath and refIndexName \n\
example: bulk.rna.seq/aale/input/ctrl.1
\n\
"


if [ ! $# -eq 3 ]; then
    printf "ERROR: usage $scriptName dataDirPath refIndexName listOfSalmonLogs \n\n"
    printf "$usage \n"
    exit 1
fi

dataDirPath=$1
ref=$2

listOfSalmonLogs=$3
if [[ ! -f $listOfSalmonLogs ]]
then
    printf "ERROR: listOfSalmonLogs $listOfSalmonLogs does not exist\n"
    exit 1
fi


# set -x # debug trace on
# set +x # debug trace off

#ref=sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
sedExpr="'s/$ref.*\$//g'"

#(base) [aedavids@mustard aale.kras]$ echo "${sedExpr}"
#'s/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.*$//g'

#eval sed "${sedExpr}" aale.kras.sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.salmon.logs.txt
eval sed "${sedExpr}" $listOfSalmonLogs | sed 's/\/private\/groups\/kimlab\/aale.kras\/data\///g'

