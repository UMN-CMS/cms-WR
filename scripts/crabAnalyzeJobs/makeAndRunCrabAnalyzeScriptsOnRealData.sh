#!/bin/bash
#all three arrays must have the same number of elements!
datasets=(DoubleMuonSkimmedMuMuJJ MuonEGSkimmedEMuJJ)
channel=(mumujj emujj)
inputData=('/DoubleMuon/skalafut-realData_DoubleMuon_13TeV_50ns_mumujj_signalAndLowMassRegionSkim_MINIAODSIM_sideband_output-8f2820a5f32d5fceefe9beeccf804df3/USER'  '/MuonEG/skalafut-realData_MuonEG_13TeV_50ns_emujj_signalAndLowMassRegionSkim_MINIAODSIM_sideband_output-8766a3daa7e3e97aac1ccb2325ba1f53/USER')

for q in ${!datasets[*]}
do

		#replace SMPL with an element from datasets, CHNL with an element from channel, and INPTDATA with an element from inputData
		#in analyze_smpl_skims_chnl_crab.py
		#inputData, datasets, and channel have the same number of elements
		eval "sed 's/SMPL/${datasets[$q]}/g' analyze_smpl_skims_chnl_crab.py > analyze_smpl_skims_chnl_crab_one.py"
		eval "sed 's/CHNL/${channel[$q]}/g' analyze_smpl_skims_chnl_crab_one.py > analyze_smpl_skims_chnl_crab_two.py"
		eval "sed 's@INPTDATA@${inputData[$q]}@g' analyze_smpl_skims_chnl_crab_two.py > analyze_${datasets[$q]}_skims_${channel[$q]}_crab.py"
		rm analyze_smpl_skims_chnl_crab_one.py analyze_smpl_skims_chnl_crab_two.py
		
		#submit jobs to crab using the newly created crab skim .py file
		#echo "crab submit -c analyze_${datasets[$q]}_skims_${channel[$q]}_crab.py"
		eval "crab submit -c analyze_${datasets[$q]}_skims_${channel[$q]}_crab.py"
	
done


