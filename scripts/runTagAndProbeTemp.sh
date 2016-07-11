#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 8nh -J DataEEtagandprobe < scripts/runDataEEtagandprobeJob.sh
UP=LOCALPATH
export X509_USER_PROXY=PROXYPATH

cd $UP
eval `scramv1 runtime -sh`
eval 'make'
eval "./bin/analysis -m MODE -c CHNL --isTagAndProbe true --ignoreDyScaleFactors UPDATE"

