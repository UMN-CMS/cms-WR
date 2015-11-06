#!/bin/bash

#lowest values for WR and Nu masses
wrMass=800
nuMass=50
minNuMass=50
maxWrMass=3500
increment=50

#number of GEN evts per job
#1000 evts at GEN lvl takes about 2 minutes
numEvts=15000
q=1

while [ $wrMass -le $maxWrMass ]; do
	while [ $nuMass -lt $wrMass ]; do
		#echo wrMass is $wrMass
		#echo nuMass is $nuMass
		#echo ""

		#replace NUM, MMAASS, and MASSNU in .csh, WR_GEN.py, and gen job submission files
		#also replace EVTCT in WR_GEN.py file

		#replacements in run...csh file
		eval "sed 's/NUM/$q/g' runGenAndOfflineAnalysis.csh > run_one.csh"
		eval "sed 's/MMAASS/$wrMass/g' run_one.csh > run_two.csh"
		eval "sed 's/MASSNU/$nuMass/g' run_two.csh > runGenAndOfflineAnalysis_WR_${wrMass}_NU_${nuMass}_${q}.csh"
		rm run_one.csh run_two.csh

		#replacements in job submission file
		eval "sed 's/NUM/$q/g' genAndRunOfflineAnalysis > job_one"
		eval "sed 's/MMAASS/$wrMass/g' job_one > job_two"
		eval "sed 's/MASSNU/$nuMass/g' job_two > genAndRunOfflineAnalysis_WR_${wrMass}_NU_${nuMass}_${q}"
		rm job_one job_two 

		#replacements in WR_GEN.py file
		eval "sed 's/NUM/$q/g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN.py > wrGen_one.py"
		eval "sed 's/MMAASS/$wrMass/g' wrGen_one.py > wrGen_two.py"
		eval "sed 's/MASSNU/$nuMass/g' wrGen_two.py > wrGen_three.py"
		eval "sed 's/EVTCT/$numEvts/g' wrGen_three.py > WR_M-${wrMass}_ToLNu_M-${nuMass}_GEN_${q}.py"
		rm wrGen_one.py wrGen_two.py wrGen_three.py


		#increment nuMass before restarting the loop
		let nuMass=nuMass+$increment
	done
	#reset the nuMass value to the lower bound value before restarting the while loop
	#with a higher wrMass value
	let nuMass=minNuMass
	let wrMass=wrMass+$increment

done

