#!/bin/bash

#use this script to run the primary python files (with cmsRun) which generate TTrees for various signal, sideband, and check regions

#analyzersToRun contains the names of the python files to be executed with cmsRun, excluding the leading test/
analyzersToRun=('checkZeeChain_cfg.py'  'check_emu_leptons.py'  'quickRecoKinematicsAnalyzer_eejj_cfg.py'  'quickRecoKinematicsAnalyzer_emujj_cfg.py')
analyzerDir='test'
stdOutDir='stdOutFromAnalysis'
inputFileDir='fileLocations'
mainOutputDir='/eos/uscms/store/user/skalafut'
workingDir='/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR'


#each of the following arrays must have the same number of elements
#NOTE the order of elements in these arrays is important
#the Nth element in one array should correspond to the Nth element in all other arrays
stdOutFileName=('analysisOutput_doubleEGCheckZee25nsRun2015D_Oct29_noHLT.txt'  'analysisOutput_muonEGCheckEMu_Run2015D25ns_Oct29_noHLT.txt'  'analysisOutput_doubleEGToEEJJLowMassSideband25nsRun2015D_Oct29_withHLT.txt'  'analysisOutput_muonEGToEMuJJLowMassSideband_Run2015D25ns_Oct29_withHLT.txt')
inputFilesList=('doubleEGHasTwoHEEPRun2015D25ns_Oct_28_formatted.txt'  'muonEGHasHEEPandIsHighPtRun2015D25ns_Oct_28_formatted.txt'  'doubleEGHasTwoHEEPRun2015D25ns_Oct_28_formatted.txt'  'muonEGHasHEEPandIsHighPtRun2015D25ns_Oct_28_formatted.txt')
outputFileDir=('analyzed_25ns_skims_check_Zee_using_two_HEEP'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj'  'analyzed_25ns_skims_low_dilepton_mass_emujj')
outputTreeFileName=('analyzed_DoubleEG_skim_25ns_hasTwoHEEP_Run2015D.root'  'analyzed_MuonEG_25ns_skim_check_emu_noHLT_Run2015D.root'  'analyzed_DoubleEG_25ns_skim_low_mass_region_eejj_Run2015D.root'  'analyzed_MuonEG_25ns_skim_low_dilepton_mass_region_emujj_Run2015D.root')
globalTag=('74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2')


for r in ${!inputFilesList[*]}
do
	#loop over all lists of input files
	for q in ${!analyzersToRun[*]}
	do
		#loop over all python analyzers/TTree producers which should be used
		#skip combinations which do not make sense, like checkZeeChain and an input list of muonEG files
		if [ ["${analyzersToRun[$q]}" == 'checkZeeChain_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_Zee_using_two_HEEP'] -o [ "${analyzersToRun[$q]}" == 'check_emu_leptons.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_emu' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_eejj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_emujj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_mass_emujj' ] ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$r]}/${outputTreeFileName[$r]}"

		fi

	done


done

