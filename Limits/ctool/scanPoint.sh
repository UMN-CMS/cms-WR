#!/bin/sh

cmswd=$1
jobwd=$2
inp=$3
mass=$4
outp=$5

method=$6
toys=$7
seed=1789

source /local/cms/sw/cmsset_default.sh

cd ${wd}

cmsenv

cd ${jobwd}

combine -v0 -t${toys} -n WRmu -m ${mass} -M ${method} -H ProfileLikelihood -s ${seed}  ${inp} > ${outp} 2>&1 

#-m 16001000 -s1002 -n WRmu -t100 -v0  