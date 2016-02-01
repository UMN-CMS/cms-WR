#!/bin/csh

#Setup environment, move and untar needed files and directories, compile everything
setenv HOME /eos/uscms/store/user/skalafut/WR/13TeV/analyzed_WRSignal_Trees
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog
setenv GENSIM /eos/uscms/store/user/skalafut/WR/13TeV/Signal_GEN-SIM

#source /uscmst1/prod/sw/cms/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh prod
setenv SCRAM_ARCH slc6_amd64_gcc491

cmsrel CMSSW_7_4_1_patch1
mv cmsswReleaseForGenSimToAnalyzedTrees.tar CMSSW_7_4_1_patch1/src/.
cd CMSSW_7_4_1_patch1/src/
tar -xvf cmsswReleaseForGenSimToAnalyzedTrees.tar
rm cmsswReleaseForGenSimToAnalyzedTrees.tar
cd ../..
echo "unzipped CMSSW_7_4_1_patch1 and cmsWR files"

mv pu.tar step2_DIGI_RECO_WR_MMAASS_NU_MASSNU_NUM.py step3_PAT_WR_MMAASS_NU_MASSNU_NUM.py WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_NUM.py CMSSW_7_4_1_patch1/src/ExoAnalysis/cmsWR/test/.

echo "moved PU files, gen-sim, and reco python files into src/ dir of CMSSW_7_4_1_patch1"

cd CMSSW_7_4_1_patch1/src/ExoAnalysis/cmsWR/test/
tar -xvf pu.tar
rm pu.tar
echo "unzipped pu tar file"

cmsenv
eval `scramv1 runtime -csh`
eval 'scram b -j 8'

echo "called cmsenv and scram b"

cmsRun WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_NUM.py
echo "finished generating evts at GEN-SIM lvl"

cmsRun step2_DIGI_RECO_WR_MMAASS_NU_MASSNU_NUM.py
mv WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_13TeV_NUM.root $GENSIM/.
echo "finished simulation up to RECO and moving GENSIM file to skalafut/WR/13TeV/Signal_GEN-SIM"

cmsRun step3_PAT_WR_MMAASS_NU_MASSNU_NUM.py
echo "finished compression to miniAOD"

cmsRun unmatched_recoElectronChannel_noCmdLineInputs_cfg.py
echo "finished running HLT and offline selection"

mv analyzed_tree_eejjSignalRegion_WR_Signal.root analyzed_tree_eejjSignalRegion_WR_MMAASS_NU_MASSNU.root
echo "renamed analyzed TTree file"

mv analyzed_tree_eejjSignalRegion_WR_MMAASS_NU_MASSNU_NUM.root $HOME/.
echo "moved analyzed TTree file to my eos space" 
