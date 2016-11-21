#!/bin/bash

WrMasses=(800 1600 2400 3200 4000 6000) 
NuMasses=(400 800 1200 1600 2000 3000)  

inputDatasetsEE=('/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-800_ToLNu-400_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-1600_ToLNu-800_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu_GEN_SIM_13TeV-2016/gnegro-WR-2400_ToLNu-1200_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER' 
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-3200_ToLNu-1600_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-4000_ToLNu-2000_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER')
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-6000_ToLNu-3000_ToEEJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER')

inputDatasetsMuMu=('/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-800_ToLNu-400_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-1600_ToLNu-800_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-2400_ToLNu-1200_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-3200_ToLNu-1600_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-4000_ToLNu-2000_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-6000_ToLNu-3000_ToMuMuJJ_RECO_13TeV-2016-9b5a96132a75d10c8d7dd85a592fb225/USER')


for h in ${!WrMasses[*]}
do

	## EEJJ
	#make python config file for miniAOD step 
	eval "sed 's/MASSWR/${WrMasses[$h]}/g' WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD.py > WR_temp1.py"                        
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp1.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_EEJJ.py"
	rm WR_temp1.py 

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD_EEJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_temp4.py"
	eval "sed 's@datasetFromDBS@${inputDatasetsEE[$h]}@g' crabConfig_WR_temp4.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_EEJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py crabConfig_WR_temp4.py


	## MuMuJJ
	#make python config file for miniAOD step 
	eval "sed 's/MASSWR/${WrMasses[$h]}/g' WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD.py > WR_temp1.py"                        
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp1.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_MuMuJJ.py"
	rm WR_temp1.py 

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD_MuMuJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_temp4.py"
 	eval "sed 's@datasetFromDBS@${inputDatasetsMuMu[$h]}@g' crabConfig_WR_temp4.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_MuMuJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py crabConfig_WR_temp4.py

done