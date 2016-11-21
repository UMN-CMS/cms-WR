#!/bin/bash

WrMasses=(800 1600 2400 3200 4000 6000) 
NuMasses=(400 800 1200 1600 2000 3000)  

inputDatasetsEE=('/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-800_ToLNu-400_ToEEJJ_GEN_SIM_13TeV-2016-2aad16ff0ab847864e355906fc4955e3/USER' 
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-1600_ToLNu-800_ToEEJJ_GEN_SIM_13TeV-2016-708ba97e6ca237382265320527720c68/USER' 
	'/WR-ToLNu_GEN_SIM_13TeV-2016/gnegro-WR-2400_ToLNu-1200_GEN_SIM_13TeV-2016-fe8dbe1d9479add77a19c2f32724a4f4/USER' 
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-3200_ToLNu-1600_ToEEJJ_GEN_SIM_13TeV-2016-3c308739612684c4a68e8f1710a50456/USER' 
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-4000_ToLNu-2000_ToEEJJ_GEN_SIM_13TeV-2016-d4a85bd98c0918ac20cd8c7328485c7c/USER' 
	'/WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016/gnegro-WR-6000_ToLNu-3000_ToEEJJ_GEN_SIM_13TeV-2016-0220955e543706145c967b7831a6ec8a/USER')

inputDatasetsMuMu=('/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-800_ToLNu-400_ToMuMuJJ_GEN_SIM_13TeV-2016-4e16d641f759f8424a52cf63661ba1bf/USER' 
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-1600_ToLNu-800_ToMuMuJJ_GEN_SIM_13TeV-2016-ad6b23c4d3b2d5be69f50c0516d9c002/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-2400_ToLNu-1200_ToMuMuJJ_GEN_SIM_13TeV-2016-3c307f2f753fe377933b0e63af85565a/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-3200_ToLNu-1600_ToMuMuJJ_GEN_SIM_13TeV-2016-573e3c036207f6558177226fffd818e3/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-4000_ToLNu-2000_ToMuMuJJ_GEN_SIM_13TeV-2016-fe00b87de76cf5b5eba9334787c26528/USER'
	'/WR-ToLNu-ToMuMuJJ_GEN_SIM_13TeV-2016/gnegro-WR-6000_ToLNu-3000_ToMuMuJJ_GEN_SIM_13TeV-2016-7ebc1ff7a2cd8e0b0666cd95fcf45513/USER')


for h in ${!WrMasses[*]}
do

	## EEJJ
	#make python config file for HLT step 
	eval "sed 's/MASSWR/${WrMasses[$h]}/g' WR_M-UNDEF_ToLNu_M-UNDEF_HLT.py > WR_temp1.py"                        
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp1.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_EEJJ.py"
	rm WR_temp1.py 

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_HLT_EEJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_temp4.py"
	eval "sed 's@datasetFromDBS@${inputDatasetsEE[$h]}@g' crabConfig_WR_temp4.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_EEJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py crabConfig_WR_temp4.py


	## MuMuJJ
	#make python config file for HLT step 
	eval "sed 's/MASSWR/${WrMasses[$h]}/g' WR_M-UNDEF_ToLNu_M-UNDEF_HLT.py > WR_temp1.py"                        
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp1.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_MuMuJJ.py"
	rm WR_temp1.py 

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_HLT_MuMuJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_temp4.py"
 	eval "sed 's@datasetFromDBS@${inputDatasetsMuMu[$h]}@g' crabConfig_WR_temp4.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_MuMuJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py crabConfig_WR_temp4.py

done
