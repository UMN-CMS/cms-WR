#!/bin/bash

WrMasses=(800 1600 2400 3200 4000 6000)
NuMasses=(400 800 1200 1600 2000 3000)


for h in ${!WrMasses[*]}
do

	eval 'crab submit -c crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_EEJJ.py'
	eval 'crab submit -c crabConfig_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN_SIM_MuMuJJ.py'

done

