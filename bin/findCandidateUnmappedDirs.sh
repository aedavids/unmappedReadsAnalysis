#!/bin/bash

#
# our scripts are written so that if output director exists
# the script will skip it
#
# some time we expect change parameter and will need to remove old output
# so we can re-run
#
# assumes that the the output from tools like salmon, STAR, fastqc are under a 'unmapped' directory

#set -x   # turn debug on
# set +x # turn debug off

dataRoot=/private/groups/kimlab/kras.ipsc/data


#
# salmon unmapped reads where split into different files
# see salmon doc for details
# unmappedType = unmapped m1, m2, m12, d
#

# print warning to stderr
>&2 printf "WARNING check output carefully it may not find everything \n" ;
#find {${dataRoot}/bulk.data,${dataRoot}/exo.data} -type d -name "*unmapped*"

for i in unmapped m1 m2 m12 d;
do

    find {${dataRoot}/bulk.data,${dataRoot}/exo.data} -type d -name $i
done

