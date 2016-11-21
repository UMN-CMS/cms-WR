#!/bin/bash

WrMasses=(800 1600 2400 3200 4000 6000)
NuMasses=(400 800 1200 1600 2000 3000)


for h in ${!WrMasses[*]}
do
	
	## EEJJ
	#make python config file for GEN-SIM step                     
	eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_EEJJ.py > WR_temp.py"   
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_EEJJ.py"
	rm WR_temp.py

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_EEJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_EEJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py


	## MuMuJJ
	#make python config file for GEN-SIM step
	eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_MuMuJJ.py > WR_temp.py"    
	eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_temp.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_MuMuJJ.py"
	rm WR_temp.py

	#make crab config file
	eval "sed 's/WR_M-UNDEF/WR_M-${WrMasses[$h]}/g' crabConfig_WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_MuMuJJ.py > crabConfig_WR_temp1.py"
	eval "sed 's/ToLNu_M-UNDEF/ToLNu_M-${NuMasses[$h]}/g' crabConfig_WR_temp1.py > crabConfig_WR_temp2.py"
	eval "sed 's/WR-MUNDEF/WR-${WrMasses[$h]}/g' crabConfig_WR_temp2.py > crabConfig_WR_temp3.py"
	eval "sed 's/ToLNu-MUNDEF/ToLNu-${NuMasses[$h]}/g' crabConfig_WR_temp3.py > crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_MuMuJJ.py"
	rm crabConfig_WR_temp1.py crabConfig_WR_temp2.py crabConfig_WR_temp3.py

done

