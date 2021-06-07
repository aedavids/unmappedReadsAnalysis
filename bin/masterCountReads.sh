#!/bin/bash
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/13/21
#

# set -x # turn debug on
# set + x # turn debug off

scriptName=$0
if [ $# -lt 1 ];
then
    printf "ERROR $scriptName: no salmonLog files passed \n"
    exit 1
fi


readFileTypes='mapped unmapped d u m1 m2 m12'

# print header
hdr=""
for t in $readFileTypes;
do
    #    printf "${t}Reads\t"
    if [ "$hdr" == "" ];
    then
        hdr="${t}"
    else
        hdr="${hdr}\t${t}"
    fi
    
done
printf "${hdr}\n"


for salmonLog in $@;
do
    countString=""
    unmappedDir=`dirname $salmonLog`/../unmapped
    if [ -d $unmappedDir ];
    then
        #printf "\nAEDWIP           unmappedDir: $unmappedDir \n"        
        for t in $readFileTypes;
        do
            if [ $t == "mapped" ] ;
            then
                mate1="${unmappedDir}/mate1.mapped.fastq"
            elif [ $t == "unmapped" ];
            then
                mate1="${unmappedDir}/mate1.unmapped.fastq"                
            else
                mate1="${unmappedDir}/mate1.${t}.unmapped.fastq"
            fi

            c=`countReads.sh ${mate1}`
            #printf "AEDWIP ${mate1} c:$c\n"            

            if [ "$countString" == "" ];
            then
                countString="${c}"
            else
                countString="${countString}\t${c}"
            fi
        done
        
    else
        #printf "\nAEDWIP not found unmappedDir: $unmappedDir \n"
        for t in $readFileTypes;
        do
            if [ "$countString" == "" ];
            then
                countString="Na"
            else
                countString="${countString}\tNa"
            fi
        done
    fi

    printf "${countString}\n"

done


