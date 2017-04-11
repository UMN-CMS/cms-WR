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
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_MuMuJJ_NUM.py"
mv WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_MuMuJJ_NUM.root /tmp/skalafut/.
rm WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_MuMuJJ_NUM.py

#run HLT on GEN-SIM root file, then delete HLT py file and GENSIM root file made in previous step
eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_HLT_MuMuJJ_NUM.py"
mv WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_MuMuJJ_NUM.root /tmp/skalafut/.
rm WR_M-MASSWR_ToLNu_M-MASSNU_HLT_MuMuJJ_NUM.py
#rm WR_M-MASSWR_ToLNu_M-MASSNU_HLT_MuMuJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_MuMuJJ_NUM.root

##run RECO on RAW-SIM root file, then delete RECO py file and HLT root file made in previous step
#eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_RECO_MuMuJJ_NUM.py"
#mv WR_M-MASSWR_ToLNu_M-MASSNU_AODSIM_13TeV-2016_MuMuJJ_NUM.root /tmp/skalafut/.
#rm WR_M-MASSWR_ToLNu_M-MASSNU_RECO_MuMuJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_MuMuJJ_NUM.root
#
##run miniAODSIM compression on AODSIM root file, then delete miniAOD py file and AODSIM root file made in previous step
#eval "cmsRun WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_MuMuJJ_NUM.py"
#mv WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_MuMuJJ_NUM.root /tmp/skalafut/.
#rm WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_MuMuJJ_NUM.py /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_AODSIM_13TeV-2016_MuMuJJ_NUM.root
#
##run WR skim on miniAODSIM file produced by previous step
#eval "cd ../../"
#eval "cmsRun test/temp_skims_cfg.py files=file:/tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_MuMuJJ_NUM.root output=WRtoMuMuJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root"
#mv WRtoMuMuJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root /tmp/skalafut/.
#rm /tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_MuMuJJ_NUM.root
#
##process WR skim file into WR minitrees
#eval "cmsRun test/temp_runAnalysis_cfg.py isMC=1 datasetTag=WRtoMuMuJJ_MASSWR_MASSNU_PrivReco files=file:/tmp/skalafut/WRtoMuMuJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root output=WRtoMuMuJJ_MASSWR_MASSNU_PrivReco_PartNUM.root"
#rm /tmp/skalafut/WRtoMuMuJJ_MASSWR_MASSNU_PrivReco_skimPartNUM.root
#
