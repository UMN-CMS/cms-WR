#!/bin/bash

#use this script to run the primary python files (with cmsRun) which generate TTrees for various signal, sideband, and check regions

#analyzersToRun contains the names of the python files to be executed with cmsRun, excluding the leading test/
#analyzersToRun=('checkZeeChain_cfg.py'  'check_emu_leptons.py'  'quickRecoKinematicsAnalyzer_eejj_cfg.py'  'quickRecoKinematicsAnalyzer_emujj_cfg.py')
analyzersToRun=('unmatched_recoElectronChannel_cfg.py')
analyzerDir='test'
stdOutDir='stdOutFromAnalysis'
inputFileDir='fileLocations'
mainOutputDir='/eos/uscms/store/user/skalafut'
workingDir='/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR'


#each of the following arrays must have the same number of elements
#NOTE the order of elements in these arrays is important
#the Nth element in one array should correspond to the Nth element in all other arrays
#the most important array is outputFileDir.  Mistakes in elements of this array will cause this script to fail.
stdOutFileName=('analysisOutput_dyJetsAMCNLOReMiniAODToEEJJSignalRegion_Nov11_M_50_with_MET_and_HLT.txt'  'analysisOutput_dyJetsMadgraphReMiniAODToEEJJSignalRegion_Nov11_M_50_with_MET_and_HLT.txt')
inputFilesList=('dyJetsAMCNLOHasTwoHEEP25nsNov11_formatted.txt'  'dyJetsMadgraphHasTwoHEEP25nsNov11_formatted.txt')
#outputFileDir=('analyzed_25ns_skims_check_Zee_using_two_HEEP'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj'  'analyzed_25ns_skims_low_dilepton_mass_emujj')
outputFileDir=('analyzed_25ns_eejj_signal_region')
outputTreeFileName=('analyzed_DYJets_amcnlo_M_50_25ns_eejj_signal_region_reMiniAOD.root'  'analyzed_DYJets_Madgraph_M_50_25ns_eejj_signal_region_reMiniAOD.root')
globalTag=('74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2')


for r in ${!inputFilesList[*]}
do
	#loop over all lists of input files
	for q in ${!analyzersToRun[*]}
	do
		#loop over all python analyzers/TTree producers which should be used
		#skip combinations which do not make sense, like checkZeeChain and an input list of muonEG files
		#if [ ["${analyzersToRun[$q]}" == 'unmatched_recoElectronChannel_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_eejj_signal_region'] -o ["${analyzersToRun[$q]}" == 'checkZeeChain_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_Zee_using_two_HEEP'] -o [ "${analyzersToRun[$q]}" == 'check_emu_leptons.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_emu' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_eejj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_emujj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_mass_emujj' ] ]; then
		if [ "${analyzersToRun[$q]}" == 'unmatched_recoElectronChannel_cfg.py' -a "${outputFileDir[$q]}" == 'analyzed_25ns_eejj_signal_region' ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$q]}/${outputTreeFileName[$r]} >& $workingDir/$stdOutDir/${stdOutFileName[$r]} &"

		fi

	done


done

