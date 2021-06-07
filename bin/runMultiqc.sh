#!/bin/bash
# aedavids@ucsc.edu
# 5/11/21



# ../../bin/runMultiqc.sh 2>&1 | tee runM.out.`~/extraCellularRNA/bin/dateStamp.sh`
# https://multiqc.info/docs/#clashing-sample-names


source createTmpFile.sh

set -x

scriptName=$0
root=/private/groups/kimlab/kras.ipsc/data

#
# --force
# overwrite existing report
#
# --dirs
# This will prefix every sample name with the directory path for that log file.
# useful for debugging, prevents samples from being overwritten
#
# -f, --force                     Overwrite any existing reports
# -d, --dirs                      Prepend directory to sample names
# -dd, --dirs-depth INTEGER       Prepend [INT] directories to sample names.
#                                 Negative number to take from start of path.

# -s, --fullnames                 Do not clean the sample names (leave as full
#                                 file name)

#        --file-list /private/home/aedavids/unmappedReadsAnalysis/data/multiqc.out/multiqc.input.lst.txt

p=`pwd`

# change directories
# path is to long, get warning 'duplicate sample names private gropus kimlab kras.ispc data bulk

#         --dirs 
#    cd ${root}

#         --fullnames 
# dirs-depth 4

# once we debug change output dir to root and code dirs searched
outputDir=${p}/${scriptName}.`~/extraCellularRNA/bin/dateStamp.sh`

#        --dirs --dirs-depth 6

#
# focus on gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
#

# ignore did not work as expeced with star and fastq, we get
# m12, m1, m2, ...
# unmappedTmp=`createTmpFile`
# add_on_exit rm  $unmappedTmp
# ./findCandidateUnmappedDirs.sh |grep -v "unmapped/d" | grep -v unmapped/m12 | grep -v unmapped/m2 | grep -v unmapped/m1 | grep unmapped/unmapped > $unmappedTmp

# head $unmappedTmp
# --file-list $unmappedTmp

# used --module to speed up. Do not check modules we are not interseted in
cd $root
multiqc ./bulk.data ./exo.data \
        --dirs \
        --fullnames \
        --verbose \
        --force \
        --module salmon \
        --module star \
        --module htseq \
        --module fastqc \
        --outdir $outputDir \
        --ignore cross.compare.de.seq/ \
        --ignore nanopore.data/ \
        --ignore single.cell.data/ \
        --ignore test.gen1c.exo/ \
        --ignore deep.variant.output/ \
        --ignore gencode.salmon.out/ \
        --ignore gencode.te.locus.salmon.out/ \
        --ignore te.combined.salmon.out \
        --ignore locus.te.combined.salmon.out/ \
        --ignore hap.py.output \
        --ignore */unmapped/d \
        --ignore */unmapped/u \
        --ignore */unmapped/m12 \
        --ignore */unmapped/m1 \
        --ignore */unmapped/m2 \
        --ignore */unmapped/mate*d*fastq \
        --ignore */unmapped/mate*u*fastq \
        --ignore */unmapped/mate*m12*fastq \
        --ignore */unmapped/mate*m1*fastq \
        --ignore */unmapped/mate*m2*fastq

