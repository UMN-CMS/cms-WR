#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 8nh -J jobName < scriptName.sh
UP=LOCALPATH
export X509_USER_PROXY=PROXYPATH

cd $UP
eval `scramv1 runtime -sh`
eval "cmsenv"
eval "cmsRun WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py"
eval "cmsRun test/checkWRDecay_crabSafe_cfg.py channel=CHNL files=file:WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_13TeV_NUM.root output=OUTPATH/analyzed_SMPL_MMAASS_NU_MASSNU_NUM.root"
rm WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_13TeV_NUM.root

