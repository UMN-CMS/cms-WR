#!/bin/bash

nuMass=(2080 520)
wrMass=(2600 2600)
dir=('/eos/uscms/store/user/skalafut/WR_2600_ToENu_2080_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-2080_13TeV_GEN-SIM/150519_203331/0000'  '/eos/uscms/store/user/skalafut/WR_2600_ToENu_520_ToEEJJ_GENSIM/WRToENuToEEJJ_MW-2600_MNu-520_13TeV_GEN-SIM/150519_203208/0000')

for x in ${!wrMass[*]}
do

	#max = 80
	for q in {1..80} 
	do #replace NUM, MMAASS, and MASSNU in .csh, step2 .py, and inputGenToAnalyzedTreeMaker job submission files
		#also replace pathToInputGenSim in inputGenToAnalyzedTreeMaker job submission files

		#replacements in run...csh file
		eval "sed 's/NUM/$q/g' runGenSimToAnalyzedTree.csh > run_one.csh"
		eval "sed 's/MMAASS/${wrMass[$x]}/g' run_one.csh > run_two.csh"
		eval "sed 's/MASSNU/${nuMass[$x]}/g' run_two.csh > runGenSimToAnalyzedTree_WR_${wrMass[$x]}_NU_${nuMass[$x]}_$q.csh"
		rm run_one.csh run_two.csh

		#replacements in step2_DIGI_RECO.py file
		eval "sed 's/NUM/$q/g' step2_DIGI_RECO.py > reco_one.py"
		eval "sed 's/MMAASS/${wrMass[$x]}/g' reco_one.py > reco_two.py"
		eval "sed 's/MASSNU/${nuMass[$x]}/g' reco_two.py > step2_DIGI_RECO_WR_${wrMass[$x]}_NU_${nuMass[$x]}_$q.py"
		rm reco_one.py reco_two.py

		#nothing in step3_PAT.py file needs to be replaced with sed
		eval "sed 's/NUM/$q/g' step3_PAT.py > step3_PAT_WR_${wrMass[$x]}_NU_${nuMass[$x]}_$q.py"

		#replacements in job submission file
		eval "sed 's/NUM/$q/g' inputGenToAnalyzedTreeMaker > job_one"
		eval "sed 's/MMAASS/${wrMass[$x]}/g' job_one > job_two"
		eval "sed 's@pathToInputGenSim@${dir[$x]}@g' job_two > job_three"
		eval "sed 's/MASSNU/${nuMass[$x]}/g' job_three > inputGenToAnalyzedTreeMaker_WR_${wrMass[$x]}_NU_${nuMass[$x]}_$q"
		rm job_one job_two job_three 

	done

done

