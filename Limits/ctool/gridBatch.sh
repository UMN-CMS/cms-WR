#!/bin/sh

cmswd=$1
jobwd=$2
inp=$3
outstub=$4

source /local/cms/sw/cmsset_CMSSW6X.sh

cd ${cmswd}

cmsenv

cd ${jobwd}

${cmswd}/gridBaseLimiter.perl ${inp} ${jobwd}/${outstub}.root ${jobwd}/${outstub}.log 

