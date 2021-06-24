#!/bin/bash

set -x


printf "\n [INFO] begin `~/extraCellularRNA/bin/dateStamp.sh` \n\n"

outdirRoot=/private/groups/kimlab/aedavids/exceRpt.out/smRNA-seqDataSet
outdir="${outdirRoot}"`~/extraCellularRNA/bin/dateStamp.sh`
mkdir -p ${outdir}
chmod a+w ${outdir}

USER_ID=`id -u`
# -v ~/DirectoryContainingMyInputSample:/exceRptInput
docker run  \
       -v ${outDir}:/exceRptOutput \
       -v /scratch/aedavids/exceRpt/db:/exceRpt_DB/hg38 \
       --user=${USER_ID}\
       -t rkitchen/excerpt \
       INPUT_FILE_PATH=/exceRptInput/testData_human.fastq.gz

printf "\n\n [INFO] end `~/extraCellularRNA/bin/dateStamp.sh` \n"
