#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 8nh -J jobName < scriptName.sh
UP=LOCALPATH
#export X509_USER_PROXY=PROXYPATH

cd $UP
eval `scramv1 runtime -sh`
eval "cmsenv"

#after each new root file is made, move the file to /tmp/skalafut to avoid problems with core dumps and exceeding local disk space limit

#make the GEN-SIM root file, then delete GEN-SIM python file
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_EEJJ_NUM.py"
#eval "cp WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/."
mv WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/.
rm WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_EEJJ_NUM.py

#run HLT on GEN-SIM root file, then delete HLT py file and GENSIM root file made in previous step
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_HLT_EEJJ_NUM.py"
#eval "cp WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/."
#rm WR_M-MASSWR_ToLNu_M-MASSNU_HLT_EEJJ_NUM.py
mv WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/.
rm WR_M-MASSWR_ToLNu_M-MASSNU_HLT_EEJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_EEJJ_NUM.root

#run RECO on RAW-SIM root file, then delete RECO py file and HLT root file made in previous step
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_RECO_EEJJ_NUM.py"
mv WR_M-MASSWR_ToLNu_M-MASSNU_AODSIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/.
#eval "cp WR_M-MASSWR_ToLNu_M-MASSNU_AODSIM_13TeV-2016_EEJJ_NUM.root /tmp/skalafut/."
rm WR_M-MASSWR_ToLNu_M-MASSNU_RECO_EEJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_EEJJ_NUM.root

#run miniAODSIM compression on AODSIM root file, then delete miniAOD py file and AODSIM root file made in previous step
#move the miniAOD .root file two directories up to facilitate skimming and minitree production on lxplus batch system
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_EEJJ_NUM.py"
mv WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_EEJJ_NUM.root ../../.
rm WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_EEJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_AODSIM_13TeV-2016_EEJJ_NUM.root

#run WR skim on miniAODSIM file produced by previous step
#dont move skim file to tmp directory, jobs running on lxplus batch system cant find /tmp/skalafut/WR_M...root
eval "cd ../../"
eval "cmsRun test/temp_skims_cfg.py files=file:WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_EEJJ_NUM.root output=WRtoEEJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root"
rm WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_EEJJ_NUM.root

#process WR skim file into WR minitrees
eval "cmsRun test/temp_runAnalysis_cfg.py test=5 isMC=1 datasetTag=WRtoEEJJ_MASSWR_MASSNU_PrivReco files=file:WRtoEEJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root output=WRtoEEJJ_MASSWR_MASSNU_PrivReco_PartNUM.root"
rm WRtoEEJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root

