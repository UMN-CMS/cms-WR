#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 8nh -J jobName < scriptName.sh
UP=LOCALPATH
#grid proxy no longer needed export X509_USER_PROXY=PROXYPATH

cd $UP
eval `scramv1 runtime -sh`
eval "cmsenv"
eval "cmsRun WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py"

#move the GEN WR root file to /tmp/skalafut to prevent job crash and core dump due to many simultaneous running jobs
eval "mv WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_8TeV_NUM.root /tmp/skalafut/."
eval "cmsRun test/checkWRDecay_crabSafe_cfg.py channel=CHNL files=file:/tmp/skalafut/WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_8TeV_NUM.root output=OUTPATH/analyzed_SMPL_MMAASS_NU_MASSNU_NUM.root"
#rm WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_8TeV_NUM.root
rm WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py /tmp/skalafut/WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_8TeV_NUM.root

