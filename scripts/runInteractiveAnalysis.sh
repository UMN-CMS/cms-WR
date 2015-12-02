#!/bin/bash

#use this script to run the primary python files (with cmsRun) which generate TTrees for various signal, sideband, and check regions

#analyzersToRun contains the names of the python files to be executed with cmsRun, excluding the leading test/
#analyzersToRun=('checkZeeChain_cfg.py'  'check_emu_leptons.py'  'quickRecoKinematicsAnalyzer_eejj_cfg.py'  'quickRecoKinematicsAnalyzer_emujj_cfg.py'  'unmatched_recoElectronChannel_cfg.py')
analyzersToRun=('check_emu_leptons.py' 'quickRecoKinematicsAnalyzer_emujj_cfg.py')
analyzerDir='test'
stdOutDir='stdOutFromAnalysis'
inputFileDir='fileLocations'
mainOutputDir='/eos/uscms/store/user/skalafut'
workingDir='/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR'


#each of the following arrays must have the same number of elements
#NOTE the order of elements in these arrays is important
#the Nth element in one array should correspond to the Nth element in all other arrays
#the most important array is outputFileDir.  Mistakes in elements of this array will cause this script to fail.
stdOutFileName=('analysisOutput_singleMuonCheckEMu_Run2015Dv3_golden_25ns_Nov24_noHLT.txt'  'analysisOutput_singleMuonCheckEMu_Run2015Dv4_golden_25ns_Nov24_noHLT.txt'  'analysisOutput_dyJetsMadgraphCheckEMu_Nov24_noHLT.txt'  'analysisOutput_ttBarCheckEMu_Nov24_noHLT.txt'  'analysisOutput_wJetsMadgraphCheckEMu_Nov24_noHLT.txt'  'analysisOutput_wzCheckEMu_Nov24_noHLT.txt'  'analysisOutput_zzCheckEMu_Nov24_noHLT.txt'  'analysisOutput_singleMuonToEMuJJLowMassSideband_Run2015Dv3_golden_25ns_Nov24_withHLT.txt'  'analysisOutput_singleMuonToEMuJJLowMassSideband_Run2015Dv4_golden_25ns_Nov24_withHLT.txt'  'analysisOutput_dyJetsMadgraphToEMuJJLowMassSideband_Nov24_withHLT.txt'  'analysisOutput_ttBarToEMuJJLowMassSideband_Nov24_withHLT.txt'  'analysisOutput_wJetsMadgraphToEMuJJLowMassSideband_Nov24_withHLT.txt'  'analysisOutput_wzToEMuJJLowMassSideband_Nov24_withHLT.txt'  'analysisOutput_zzToEMuJJLowMassSideband_Nov24_withHLT.txt')
inputFilesList=('singleMuonHasHEEPandIsHighPtRun2015Dv3_25ns_gold_Nov_24_formatted.txt'  'singleMuonHasHEEPandIsHighPtRun2015Dv4_25ns_gold_Nov_24_formatted.txt'  'dyJetsMadgraphHasHEEPandIsHighPtID25nsNov11_formatted.txt'  'ttBarPowhegPythiaReMiniAODToEMuHasHEEPandIsHighPt25ns_Oct16_formatted.txt'  'wJetsMadgraphReMiniAODHasHEEPandIsHighPtID25ns_Nov24_formatted.txt'  'wzHasHEEPandIsHighPtID25ns_Nov24_formatted.txt'  'zzHasHEEPandIsHighPtID25ns_Nov24_formatted.txt'  'singleMuonHasHEEPandIsHighPtRun2015Dv3_25ns_gold_Nov_24_formatted.txt'  'singleMuonHasHEEPandIsHighPtRun2015Dv4_25ns_gold_Nov_24_formatted.txt'  'dyJetsMadgraphHasHEEPandIsHighPtID25nsNov11_formatted.txt'  'ttBarPowhegPythiaReMiniAODToEMuHasHEEPandIsHighPt25ns_Oct16_formatted.txt'  'wJetsMadgraphReMiniAODHasHEEPandIsHighPtID25ns_Nov24_formatted.txt'  'wzHasHEEPandIsHighPtID25ns_Nov24_formatted.txt'  'zzHasHEEPandIsHighPtID25ns_Nov24_formatted.txt')

###SAVE THIS list of all output file dirs outputFileDir=('analyzed_25ns_skims_check_Zee_using_two_HEEP'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj'  'analyzed_25ns_skims_low_dilepton_mass_emujj')

outputFileDir=('analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_check_emu'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj'  'analyzed_25ns_skims_low_dilepton_mass_emujj')
outputTreeFileName=('analyzed_SingleMuon_25ns_skim_check_emu_noHLT_Run2015D_v3_golden.root'  'analyzed_SingleMuon_25ns_skim_check_emu_noHLT_Run2015D_v4_golden.root'  'analyzed_DYJets_Madgraph_25ns_skim_check_emu_noHLT.root'  'analyzed_TTOnly_PowhegPythia_25ns_skim_check_emu_noHLT_reMiniAOD.root'  'analyzed_WJets_25ns_skim_check_emu_noHLT.root'  'analyzed_WZ_25ns_skim_check_emu_noHLT.root'  'analyzed_ZZ_25ns_skim_check_emu_noHLT.root'  'analyzed_SingleMuon_25ns_skim_low_dilepton_mass_region_emujj_Run2015D_v3_golden.root'  'analyzed_SingleMuon_25ns_skim_low_dilepton_mass_region_emujj_Run2015D_v4_golden.root'  'analyzed_DYJets_Madgraph_25ns_skim_low_dilepton_mass_region_emujj_reMiniAOD.root'  'analyzed_TTOnly_PowhegPythia_25ns_skim_low_dilepton_mass_region_emujj_reMiniAOD.root'  'analyzed_WJets_25ns_skim_low_dilepton_mass_region_emujj.root'  'analyzed_WZ_25ns_skim_low_dilepton_mass_region_emujj.root'  'analyzed_ZZ_25ns_skim_low_dilepton_mass_region_emujj.root')
globalTag=('74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2'  '74X_mcRun2_asymptotic_v2')

#stdOutFileName=('analysisOutput_doubleEGToEEJJLowMassSideband2015C_Nov16_HLT.txt'  'analysisOutput_muonEGToEMuJJLowMassSideband2015C_Nov16_HLT.txt')
#inputFilesList=('doubleEGHasTwoHEEPRun2015C25ns_Oct_23_formatted.txt'  'muonEGHasHEEPandIsHighPtRun2015C25ns_Oct_23_formatted.txt')
#outputFileDir=('analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj'  'analyzed_25ns_skims_low_dilepton_mass_emujj')
#outputTreeFileName=('analyzed_DoubleEG_25ns_skim_low_mass_region_eejj_Run2015C.root'  'analyzed_MuonEG_25ns_skim_low_dilepton_mass_region_emujj_Run2015C.root')
#globalTag=('74X_dataRun2_Prompt_v2'  '74X_dataRun2_Prompt_v2')


for r in ${!inputFilesList[*]}
do
	#loop over all lists of input files
	for q in ${!analyzersToRun[*]}
	do
		#loop over all python analyzers/TTree producers which should be used
		#skip combinations which do not make sense, like checkZeeChain and an input list of muonEG files
		#if [ ["${analyzersToRun[$q]}" == 'unmatched_recoElectronChannel_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_eejj_signal_region'] -o ["${analyzersToRun[$q]}" == 'checkZeeChain_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_Zee_using_two_HEEP'] -o [ "${analyzersToRun[$q]}" == 'check_emu_leptons.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_emu' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_eejj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj' ] -o [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_emujj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_mass_emujj' ] ]; then
		if [ "${analyzersToRun[$q]}" == 'unmatched_recoElectronChannel_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_eejj_signal_region' ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$r]}/${outputTreeFileName[$r]} >& $workingDir/$stdOutDir/${stdOutFileName[$r]} &"

		fi

		if [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_eejj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_and_fourObj_mass_eejj' ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$r]}/${outputTreeFileName[$r]} >& $workingDir/$stdOutDir/${stdOutFileName[$r]} &"

		fi

		if [ "${analyzersToRun[$q]}" == 'quickRecoKinematicsAnalyzer_emujj_cfg.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_low_dilepton_mass_emujj' ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$r]}/${outputTreeFileName[$r]} >& $workingDir/$stdOutDir/${stdOutFileName[$r]} &"

		fi

		if [ "${analyzersToRun[$q]}" == 'check_emu_leptons.py' -a "${outputFileDir[$r]}" == 'analyzed_25ns_skims_check_emu' ]; then
			#call cmsRun with the appropriate cmd line options
			eval "cmsRun $workingDir/$analyzerDir/${analyzersToRun[$q]} GT=${globalTag[$r]} files=$(cat $workingDir/$inputFileDir/${inputFilesList[$r]}) output=$mainOutputDir/${outputFileDir[$r]}/${outputTreeFileName[$r]} >& $workingDir/$stdOutDir/${stdOutFileName[$r]} &"

		fi



	done


done

