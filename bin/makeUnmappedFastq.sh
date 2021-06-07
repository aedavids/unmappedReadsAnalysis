#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/7/21
# assume salmon has been run with --writeUnmappedNames
# splits up the original paired read fastqs into the following files
#
# all unmapped
# all mapped
# unmapped reads split by flag
# d   = mapping type was determined as mapping to a decoy sequence
# u   = The entire pair was unmapped. No mappings were found for either the left or right read.
# m1  = Left orphan (mappings were found for the left (i.e. first) read, but not the right).
# m2  = Right orphan (mappinds were found for the right read, but not the left).
# m12 = Left and right orphans. Both the left and right read mapped, but never to the same transcript.

set -x # turn debug on
# set +x # turn debug off

source createTmpFile.sh

# send a txt msg to script finishes
extraBin=/private/home/aedavids/extraCellularRNA/bin
add_on_exit "${extraBin}/dataIsUpSMS.sh 6508662639@txt.att.net"

# redirect file permission errors to /dev/null
listOfSalmonOutDir=`find /private/groups/kimlab/kras.ipsc/ -type d -name gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx -print 2> /dev/null`

# create an variable of type in
declare -i i
i=0

function getUncompressedFastQ() {
    #
    # if we have to make multiple passes over a fastq and are running
    # SDD it is probably faster to uncompress 1 and read from disk each time
    #
    mate=$1
    file ${mate} | grep gzip > /dev/null 
    grepExitStatus=$? # zero means found
    isGZIP=0
    if [ $grepExitStatus -eq 0 ]; then
        mateTmp=`createTmpFile`

        # can not call add_on_exit from a function.
        # it does not work
        # add_on_exit rm ${mateTmp}
        
        zcat ${mate} > ${mateTmp}
        echo $mateTmp # bash return value

        isGZIP=1
    else
        echo $mate # bash return value
    fi

    # caller can use exit status to know if
    # temp file was created and will need to be cleaned up
    return $isGZIP
}

#for salmonOut in '/private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx';
for salmonOut in $listOfSalmonOutDir;
do
    output="${salmonOut}/unmapped"
    if [ -d $output ]; then
        printf "WARN skipping $output. it already exists\n"

    else
        printf "\n\n##############\n"
        mkdir -p $output
        
        unmappedNames="${salmonOut}/aux_info/unmapped_names.txt"

        # find paired fastq file
        cmd_info="${salmonOut}/cmd_info.json"
        mates1=`grep mates1 ${cmd_info} | cut -d ':' -f 2 | tr -d "\"," `
        mates2=`grep mates2 ${cmd_info} | cut -d ':' -f 2 | tr -d "\"," `

        #
        # working with uncompressed files will save a lot of time
        #
        printf "before mates1:xxx${mates1}xxx \n"
        file ${mates1}
        printf "before mates2:xxx${mates2}xxx \n"
        file ${mates2}
        mates1=$(getUncompressedFastQ $mates1)
        isM1Tmp=$?
        if [ $isM1Tmp -ne 0 ];
        then
            add_on_exit rm ${mates1}
        fi
        
            
        mates2=$(getUncompressedFastQ $mates2)
        isM2Tmp=$?
        if [ $isM2Tmp -ne 0 ];
        then
            add_on_exit rm ${mates2}
        fi
        
        printf "after AEDWIP mates1: $mates1 mates2: $mates2 \n"
        
        # select all the unmapped reads
        # unmapped.unmapped is a strange name how ever it makes
        # down stream analysis easier because all the fastq files
        # use the same naming convension
        mate1Unmapped="${output}/mate1.unmapped.unmapped.fastq"
        selectUnmappedReadsFromFastq.sh $mates1 $unmappedNames > "${mate1Unmapped}"
        mate2Unmapped="${output}/mate2.unmapped.unmapped.fastq"        
        selectUnmappedReadsFromFastq.sh $mates2 $unmappedNames > "${mate2Unmapped}"

        #
        # split the unmapped reads fastq files by flag type
        # run these concurrently in the background
        #
        validFlags='d u m1 m2 m12'
        
        for f in $validFlags;
        do
            i=0
            for mate in {$mate1Unmapped,$mate2Unmapped};
            do
                i=${i}+1
                uOut="${output}/mate${i}.${f}.unmapped.fastq"

                #
                # select flag  runs on the smaller unmapped fastq file
                # run in background
                #
                selectFlaggedUnmappedReadsFromFastq.sh $mate $f $unmappedNames > ${uOut} &
            done
            
        done

        #
        # select all the mapped reads
        # this is slow because the fastq file is big
        # we use this to throttle the script so that other users
        # can still make progress
        #
        selectMappedReadsFromFastq.sh $mates1 $unmappedNames > "${output}/mate1.mapped.fastq"
        selectMappedReadsFromFastq.sh $mates2 $unmappedNames > "${output}/mate2.mapped.fastq"

        # clean up, so we do not fill up /tmp
        if [ $isM1Tmp -ne 0 ];
        then
            'rm' ${mates1}
        fi
        
        if [ $isM2Tmp -ne 0 ];
        then
            'rm' ${mates2}
        fi
        
    fi
    
done
