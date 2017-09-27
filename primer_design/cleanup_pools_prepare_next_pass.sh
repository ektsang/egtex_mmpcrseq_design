#!/bin/bash

# author: Emily Tsang

set -o nounset -o pipefail

# take the pools designed in a previous pass and consolodate them
# only keep pools with at least 6 sites
# for the rest, recuperate the sites and include them with the sites for which not design was possible

## The variable in caps on the RHS should be set globally, e.g. in bashrc ##
dir=${EGTEXDIR}/primerDesign
#################################################################

if [ $# -ne 2 ]; then
    echo "usage: cleanup_pools_prepare_next_pass.sh old_pass_nam new_pass_name"
    exit
fi

pass=$1
newpass=$2

poolsOutput=${dir}/${pass}.pools.txt
sitesOutput=${dir}/selected.variants.${newpass}.master.txt

# delete ouput files if they already exist
if [ -e $poolsOutput ]; then
    rm $poolsOutput
fi
if [ -e $sitesOutput ]; then
    rm $sitesOutput
fi

for passdir in ${dir}/${pass}*
do
    echo $passdir
    poolFile=${passdir}/pools.txt
    # deal with pools to keep
    keep=`cut -f1 $poolFile | uniq -c | awk '$1 >= 6 {print $2}'`
    echo $keep | sed 's/ /\n/g' | grep -F -f - $poolFile >> $poolsOutput
    # deal with sites from pools that were too small
    # (identify them and get then from that pass's site file)
    echo $keep | sed 's/ /\n/g' | grep -v -F -f - $poolFile | \
	cut -f2 | grep -F -f - ${passdir}/selected.variants.*.txt >> $sitesOutput
    # deal with sites for which no design was possible
    cat ${passdir}/remaining.sites.txt >> $sitesOutput
done

# shuffle sites output file
sort -R $sitesOutput > tmp
mv tmp $sitesOutput
