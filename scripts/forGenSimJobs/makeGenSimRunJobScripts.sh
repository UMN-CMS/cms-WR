#!/bin/bash

nuMass=1000
wrMass=2000

for q in {1..200} 
do

	eval "sed 's/NUM/$q/g' runGenSimMaker.csh > runGenSimMaker_temp_1.csh"
	eval "sed 's/MASSNU/$nuMass/g' runGenSimMaker_temp_1.csh > runGenSimMaker_1_$q.csh"
	eval "sed 's/MMAASS/$wrMass/g' runGenSimMaker_1_$q.csh > runGenSimMaker_$q.csh"
	
	eval "sed 's/NUM/$q/g' genSimMaker > genSimMaker_temp"
	eval "sed 's/MASSNU/$nuMass/g' genSimMaker_temp > genSimMaker_1_temp"
	eval "sed 's/MMAASS/$wrMass/g' genSimMaker_1_temp > genSimMaker_$q"

	eval "sed 's/NUM/$q/g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM.py > WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_$q.py"
	eval "sed 's/MASSNU/$nuMass/g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_$q.py > WR_M-UNDEF_ToLNu_M-${nuMass}_GEN_SIM_$q.py"
	eval "sed 's/MMAASS/$wrMass/g' WR_M-UNDEF_ToLNu_M-${nuMass}_GEN_SIM_$q.py > WR_M-${wrMass}_ToLNu_M-${nuMass}_GEN_SIM_$q.py"

	rm runGenSimMaker_temp_1.csh runGenSimMaker_1_$q.csh
	rm genSimMaker_temp genSimMaker_1_temp
	rm WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_$q.py WR_M-UNDEF_ToLNu_M-${nuMass}_GEN_SIM_$q.py

done


