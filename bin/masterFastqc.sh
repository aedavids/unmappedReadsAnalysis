#!/bin/bash
# aedavids@ucsc.edu
# 5/17/21

#
# run fastq on all the unmapped files
#

set -x
for unmappedType in d u m1 m2 m12 unmapped;
do
    runSalmonUnmappedFastqc.sh $unmappedType > runSalmonUnmappedFastqc.sh.$unmappedType.out.`~/extraCellularRNA/bin/dateStamp.sh` 2>&1 &
    runSTARUnmappedFastqc.sh $unmappedType   > runSTARUnmappedFastqc.sh.$unmappedType.out.`~/extraCellularRNA/bin/dateStamp.sh` 2>&1 &
done
