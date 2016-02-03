#!/bin/csh

#Setup environment, move and untar needed files and directories, compile everything
setenv HOME TREEDIR
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog
setenv GENEVTS GENDIR

source /cvmfs/cms.cern.ch/cmsset_default.csh prod
setenv SCRAM_ARCH slc6_amd64_gcc491

cmsrel CMSSW_7_4_12_patch4
mv offlineWrAnalysisCode.tar CMSSW_7_4_12_patch4/src/.
cd CMSSW_7_4_12_patch4/src/
tar -xvf offlineWrAnalysisCode.tar
rm offlineWrAnalysisCode.tar
cd ../..
echo "unzipped cmsWR files"

mv WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR/test/.
echo "moved gen python file into src/ExoAnalysis/cmsWR/test/ dir of CMSSW_7_4_12_patch4"

cd CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR/

cmsenv
eval `scramv1 runtime -csh`
eval 'scram b -j 8'
cd test/
echo "called cmsenv and scram b"

cmsRun WR_M-MMAASS_ToLNu_M-MASSNU_GEN_NUM.py
cp WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_13TeV_NUM.root $GENEVTS/.
echo "finished generating evts at GEN lvl and copying evts to disk"

#now rename the WR_GEN root file to interface with checkWRDecay_crabSafe_cfg.py
mv WR_MWR_MMAASS_ToLNu_MNu_MASSNU_GEN_13TeV_NUM.root WR_GEN_13TeV.root

cmsRun checkWRDecay_crabSafe_cfg.py
echo "finished running offline selection"

mv analyzed_genWrToEEJJFullOfflineAnalysis.root analyzed_genWrToEEJJFullOfflineAnalysis_WR_MMAASS_NU_MASSNU_NUM.root
echo "renamed analyzed TTree file"

mv analyzed_genWrToEEJJFullOfflineAnalysis_WR_MMAASS_NU_MASSNU_NUM.root $HOME/.
echo "moved analyzed TTree file to my eos space" 
