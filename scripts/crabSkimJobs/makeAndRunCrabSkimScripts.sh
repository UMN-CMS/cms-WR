#!/bin/bash
datasets=(MuonEG DoubleEG SingleMuon SingleElectron DoubleMuon)
#datasets=(MuonEG)
channel=(emujj eejj mumujj)


for q in ${!datasets[*]}
do
	
	for h in ${!channel[*]}
	do

		#replace FNLST with an element from datasets, CHNL with an element from channel
		#in skim_realData_chnl_temp.py
		eval "sed 's/FNLST/${datasets[$q]}/g' skim_realData_chnl_temp.py > skim_realData_chnl_one.py"
		eval "sed 's/CHNL/${channel[$h]}/g' skim_realData_chnl_one.py > skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py"
		rm skim_realData_chnl_one.py
		
		#these if statements delete crab skim cfg files for skims which should not be run
		#like the eejj skim over the MuonEG dataset
		if [ "${datasets[$q]}" == 'MuonEG' -a "${channel[$h]}" == 'eejj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi
		if [ "${datasets[$q]}" == 'MuonEG' -a "${channel[$h]}" == 'mumujj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi

		if [ "${datasets[$q]}" == 'SingleMuon' -a "${channel[$h]}" == 'eejj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi

		if [ "${datasets[$q]}" == 'SingleElectron' -a "${channel[$h]}" == 'mumujj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi

		if [ "${datasets[$q]}" == 'DoubleMuon' -a "${channel[$h]}" == 'eejj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi
		if [ "${datasets[$q]}" == 'DoubleMuon' -a "${channel[$h]}" == 'emujj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi

		if [ "${datasets[$q]}" == 'DoubleEG' -a "${channel[$h]}" == 'mumujj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi
		if [ "${datasets[$q]}" == 'DoubleEG' -a "${channel[$h]}" == 'emujj' ]; then
			rm skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py
		fi

		#submit jobs to crab using the newly created crab skim .py file
		#echo "crab submit -c skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py"
		eval "crab submit -c skim_realData_signalAndLowMassRegions_${channel[$h]}_${datasets[$q]}.py"


	done

done


