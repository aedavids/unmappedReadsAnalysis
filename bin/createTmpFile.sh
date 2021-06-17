#!/bin/bash
#
# aedavids@ucsc.edu
# 4/30/21
# ref: https://www.linuxjournal.com/content/use-bash-trap-statement-cleanup-temporary-files
#
# this other trick works how ever you have to keep track of the file discriptors to the
# temp file by number
# and redirect output into the file handler by number
# https://unix.stackexchange.com/questions/181937/how-create-a-temporary-file-in-shell-script
#

# list of public fuctions
# debug_on_exit
# add_on_exit
#   see example bellow
# createTmpFile

# -a type is array
declare -a on_exit_items


function on_exit {
    for i in "${on_exit_items[@]}"
    do
        #echo "!!!!!! PID: $$ AEDWIP on_exit(): $i" 1>&2 2>errors.txt
        eval $i
    done
}

function debug_on_exit {
    for i in "${on_exit_items[@]}"
    do
        echo "!!!!!! debug on_exit(): $i"
    done
    
}

#
# add_on_exit can not be called from a function. it must be called from the
# main body of script
# Example of how to create a function that that uses temp files
#
# source createTmpFile.sh
# function getUncompressedFastQ() {
#     mate=$1
#     file ${mate} | grep gzip > /dev/null 
#     grepExitStatus=$? # zero means found
#     isGZIP=0
#     if [ $grepExitStatus -eq 0 ]; then
#         mateTmp=`createTmpFile`

#         # can not call add_on_exit from a function.
#         # it does not work
#         # add_on_exit rm ${mateTmp}
        
#         zcat ${mate} > ${mateTmp}
#         echo $mateTmp # bash return value

#         isGZIP=1
#     else
#         echo $mate # bash return value
#     fi

#     # caller can use exit status to know if
#     # temp file was created and will need to be cleaned up
#     return $isGZIP
# }
#
# example of how to call getUncompressed
#
# mates2=$(getUncompressedFastQ $mates2)
# isM2Tmp=$?
# if [ $isM2Tmp -ne 0 ];
# then
#     add_on_exit rm ${mates2}
# fi

function add_on_exit {
    local n=${#on_exit_items[*]}
    on_exit_items[$n]="$*"
    if [[ $n -eq 0 ]]; then
        #echo "Setting trap" 1>&2 2>errors.txt
        trap on_exit EXIT
    fi
}



function createTmpFile  {
    #set -x
    # {} executes in process
    # () executes as sub process
    # create a tmp output file and make sure it will be
    # deleted
    #https://unix.stackexchange.com/a/181938
    prefix=`basename $0`
    # gi.ucsc.edu best practics is to use /data/tmp
    # /tmp is small
    tmp=$(mktemp /data/tmp/${prefix}.XXXXXX)
    
    # return name of tmp file
    echo $tmp
}
