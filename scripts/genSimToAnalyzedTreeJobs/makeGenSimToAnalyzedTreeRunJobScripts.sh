#!/bin/bash

nuMass=(500 900 1700 2100 2500)
wrMass=(2600)
#maxCount=300

#number of GENSIM evts per job
#each evt at GENSIM lvl takes about 1 minute
numEvts=100

for r in ${!wrMass[*]}
do

	for w in ${!nuMass[*]}
	do

		for q in {1..300} 
		do
			#replace NUM, MMAASS, and MASSNU in .csh, step2 and step3 .py, and WR_GEN_SIM.py, and genSimToAnalyzedTreeMaker job submission files
			#also replace EVTCT in WR_GEN_SIM.py file

			#replacements in run...csh file
			eval "sed 's/NUM/$q/g' runGenSimToAnalyzedTree.csh > run_one.csh"
			eval "sed 's/MMAASS/${wrMass[$r]}/g' run_one.csh > run_two.csh"
			eval "sed 's/MASSNU/${nuMass[$w]}/g' run_two.csh > runGenSimToAnalyzedTree_WR_${wrMass[$r]}_NU_${nuMass[$w]}_$q.csh"
			rm run_one.csh run_two.csh

			#replacements in step2_DIGI_RECO.py file
			eval "sed 's/NUM/$q/g' step2_DIGI_RECO.py > reco_one.py"
			eval "sed 's/MMAASS/${wrMass[$r]}/g' reco_one.py > reco_two.py"
			eval "sed 's/MASSNU/${nuMass[$w]}/g' reco_two.py > step2_DIGI_RECO_WR_${wrMass[$r]}_NU_${nuMass[$w]}_$q.py"
			rm reco_one.py reco_two.py

			#replacements in step3_PAT.py file
			eval "sed 's/NUM/$q/g' step3_PAT.py > miniAOD_one.py"
			eval "sed 's/MMAASS/${wrMass[$r]}/g' miniAOD_one.py > miniAOD_two.py"
			eval "sed 's/MASSNU/${nuMass[$w]}/g' miniAOD_two.py > step3_PAT_WR_${wrMass[$r]}_NU_${nuMass[$w]}_$q.py"
			rm miniAOD_one.py miniAOD_two.py

			#replacements in job submission file
			eval "sed 's/NUM/$q/g' genSimToAnalyzedTreeMaker > job_one"
			eval "sed 's/MMAASS/${wrMass[$r]}/g' job_one > job_two"
			eval "sed 's/MASSNU/${nuMass[$w]}/g' job_two > genSimToAnalyzedTreeMaker_WR_${wrMass[$r]}_NU_${nuMass[$w]}_$q"
			rm job_one job_two 

			#replacements in WR_GEN_SIM.py file
			eval "sed 's/NUM/$q/g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM.py > wrGenSim_one.py"
			eval "sed 's/MMAASS/${wrMass[$r]}/g' wrGenSim_one.py > wrGenSim_two.py"
			eval "sed 's/MASSNU/${nuMass[$w]}/g' wrGenSim_two.py > wrGenSim_three.py"
			eval "sed 's/EVTCT/$numEvts/g' wrGenSim_three.py > WR_M-${wrMass[$r]}_ToLNu_M-${nuMass[$w]}_GEN_SIM_$q.py"
			rm wrGenSim_one.py wrGenSim_two.py wrGenSim_three.py

		done

	done

done

