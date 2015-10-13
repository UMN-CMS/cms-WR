#!/bin/bash

#all arrays must have the same number of elements
#datasets=(SingleMuon  SingleMuon  SingleMuon  DoubleEG  DoubleEG  DoubleEG  MuonEG  MuonEG  MuonEG)
#identifier=(Run2015B-17Jul2015-v1  Run2015B-05Aug2015-v1  Run2015C-PromptReco-v1  Run2015B-17Jul2015-v1  Run2015B-05Aug2015-v1  Run2015C-PromptReco-v1  Run2015B-17Jul2015-v1  Run2015B-05Aug2015-v1  Run2015C-PromptReco-v1)
#channel=(emujj  emujj  emujj  eejj  eejj  eejj  emujj  emujj  emujj)
#suffix=('OneHEEPOneIsHighPtID_17Jul2015_50ns_Sept8_2015' 'OneHEEPOneIsHighPtID_5Aug2015_50ns_Sept8_2015' 'OneHEEPOneIsHighPtID_25ns_Sept8_2015' 'TwoHEEPIDEles_17Jul2015_50ns_Sept8_2015' 'TwoHEEPIDEles_5Aug2015_50ns_Sept8_2015' 'TwoHEEPIDEles_25ns_Sept8_2015' 'OneHEEPOneIsHighPtID_17Jul2015_50ns_Sept8_2015' 'OneHEEPOneIsHighPtID_5Aug2015_50ns_Sept8_2015' 'OneHEEPOneIsHighPtID_25ns_Sept8_2015')
#lumiFiles=('/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt')

datasets=(SingleMuon  DoubleEG  MuonEG)
identifier=(Run2015B-05Aug2015-v1  Run2015B-05Aug2015-v1  Run2015B-05Aug2015-v1)
channel=(emujj  eejj  emujj)
suffix=('OneHEEPOneIsHighPtID_5Aug2015_50ns_CorrectGTandRelease_Sept8_2015' 'TwoHEEPIDEles_5Aug2015_50ns_CorrectGTandRelease_Sept8_2015' 'OneHEEPOneIsHighPtID_5Aug2015_50ns_CorrectGTandRelease_Sept8_2015')
lumiFiles=('/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt')


for q in ${!datasets[*]}
do

		#replace FNLST with an element from datasets, CHNL with an element from channel,
		#and TAG with an element from identifier in skim_realData_chnl_temp.py
		eval "sed 's/FNLST/${datasets[$q]}/g' skim_realData_chnl_temp.py > skim_realData_chnl_one.py"
		eval "sed 's@TAG@${identifier[$q]}@g' skim_realData_chnl_one.py > skim_realData_chnl_two.py"
		eval "sed 's@UNIQUE@${suffix[$q]}@g' skim_realData_chnl_two.py > skim_realData_chnl_three.py"
		eval "sed 's@LUMI@${lumiFiles[$q]}@g' skim_realData_chnl_three.py > skim_realData_chnl_four.py"
		eval "sed 's/CHNL/${channel[$q]}/g' skim_realData_chnl_four.py > skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"
		
		rm skim_realData_chnl_one.py skim_realData_chnl_two.py skim_realData_chnl_three.py skim_realData_chnl_four.py
		
		#these if statements delete crab skim cfg files for skims which should not be run
		#this is saved just as a reference example
		#if [ "${datasets[$q]}" == 'MuonEG' -a "${channel[$q]}" == 'eejj' ]; then
		#	rm skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py
		#fi

		#submit jobs to crab using the newly created crab skim .py file
		#echo "crab submit -c skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"
		eval "crab submit -c skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"

done


