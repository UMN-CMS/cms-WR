#!/bin/bash
#use this file with a command like
#bsub -R "pool>1500" -q 1nh -J uniqueJobName < runMassWindowScan_MWR_WRMASS_CHNL.sh

cd WORKDIR
eval `scramv1 runtime -sh`
eval "python scripts/limit_window_scan_MWR_WRMASS_CHNL.py > limitResultsFinerBinningTTMCWithoutSyst_MWR_WRMASS_CHNL.txt"

