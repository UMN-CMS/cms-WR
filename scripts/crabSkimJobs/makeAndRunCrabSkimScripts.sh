#!/bin/bash

#all arrays must have the same number of elements
datasets=(DoubleEG  DoubleEG)
identifier=(Run2015D-PromptReco-v4  Run2015D-PromptReco-v3)
channel=(eejj  eejj)
suffix=('OneHEEPIDEleAndLooseOrTightDoubleEleHLT_Run2015D_v4_25ns_GoldenJSON_Nov12_2015'  'OneHEEPIDEleAndLooseOrTightDoubleEleHLT_Run2015D_v3_25ns_GoldenJSON_Nov12_2015')
lumiFiles=('/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt')

#datasets=(DoubleEG  MuonEG)
#identifier=(Run2015D-PromptReco-v3  Run2015D-PromptReco-v3)
#channel=(eejj  emujj)
#suffix=('TwoHEEPIDEles_Run2015D_v3_25ns_ExclusiveSilverJSON_Nov04_2015'  'OneHEEPOneIsHighPtID_Run2015D_v3_25ns_ExclusiveSilverJSON_Nov04_2015')
#lumiFiles=('/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_256729-258313_13TeV_PromptReco_Collisions15_25ns_JSON_ExclusiveSilver.txt' '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_256729-258313_13TeV_PromptReco_Collisions15_25ns_JSON_ExclusiveSilver.txt')



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


