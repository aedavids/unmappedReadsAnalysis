#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 6/22/21
#
# create sample names for use with mineSalmonLogs.sh

usage="usage $scriptName listOfSalmonLogs  \n\
creates a list of sample names for use with  mineSalmonLogs.sh \n\
\n
"

scriptName=`basename $0`
if [ ! $# -eq 1 ]; then
    printf "ERROR: usage $scriptName listOfSalmonLogs \n\n"
    printf "$usage \n"
    exit 1
fi

listOfSalmonLogs=$1
if [[ ! -f $listOfSalmonLogs ]]
then
    printf "ERROR: listOfSalmonLogs $listOfSalmonLogs does not exist\n"
    exit 1
fi


# set -x # debug trace on
# set +x # debug trace off

ref=sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
sedExpr="'s/$ref.*\$//g'"

#(base) [aedavids@mustard aale.kras]$ echo "${sedExpr}"
#'s/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.*$//g'

#eval sed "${sedExpr}" aale.kras.sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.salmon.logs.txt
eval sed "${sedExpr}" $listOfSalmonLogs | sed 's/\/private\/groups\/kimlab\/aale.kras\///g'

