#!/bin/bash
#all three arrays must have the same number of elements!
#datasets=(TTBarSkimmedEEJJ TTBarSkimmedEMuJJ TTBarSkimmedMuMuJJ DYJetsSkimmedEEJJ DYJetsSkimmedMuMuJJ)
#channel=(eejj emujj mumujj eejj mumujj)
#inputData=('/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-TTJets_13TeV_25ns_skim_low_mass_region_EEJJ-4bd5cd38343a61b681446d6a3d70e0de/USER'  '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-TTJets_13TeV_25ns_skim_low_mass_and_signal_regions_EMuJJ_MINIAODSIM_sideband_output-8766a3daa7e3e97aac1ccb2325ba1f53/USER'  '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-TTJets_13TeV_25ns_skim_low_mass_region_MuMuJJ-9612eb022144ded38f95abdea6254832/USER'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-DYJets_M_50_largeDataset_13TeV_25ns_skim_low_mass_region_EEJJ-4bd5cd38343a61b681446d6a3d70e0de/USER'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-DYJets_M_50_largeDataset_13TeV_25ns_skim_low_mass_region_MuMuJJ-9612eb022144ded38f95abdea6254832/USER')
datasets=(TTBarSkimmedEEJJ)
channel=(eejj)
inputData=('/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skalafut-TTJets_13TeV_25ns_skim_low_mass_region_EEJJ-4bd5cd38343a61b681446d6a3d70e0de/USER')

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
		echo "crab submit -c analyze_${datasets[$q]}_skims_${channel[$q]}_crab.py"
		#eval "crab submit -c analyze_${datasets[$q]}_skims_${channel[$q]}_crab.py"
	
done


