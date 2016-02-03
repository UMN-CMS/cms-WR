#!/bin/bash

#lowest values for WR and Nu masses
#minNuMass and nuMass should be equal
#nuMass is used as an iterator which is frequently changed, while minNuMass is used
#as a reset point for nuMass when wrMass is incremented
#the value of wrMass set here determines the lowest WR mass which will be simulated
wrMass=800
nuMass=50
minNuMass=50
maxWrMass=3500
increment=50

#number of GEN evts per job
#1000 evts at GEN lvl takes about 2 minutes
#q is a number appended to the end of the output .root file names, and has no effect on
#the event generation or subsequent offline GEN analysis
numEvts=15000
q=1

#absolute path names to output file directories
#where GEN events and TTree files will be stored
pathToTrees='/uscms/home/skalafut/nobackup/WR_starting2015/genWrTTrees'
pathToGenEvts='/eos/uscms/store/user/skalafut/WR/13TeV/WRSignal_slimmedGEN'

#absolute path name to the directory where the job executable .csh file is stored
pathToExe='/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR/scripts/genAndFullOfflineAnalysisJobs'

####
##NOTE these 6 sed commands are added to maintain functionality of the root macro
##which calculates the GEN efficiency values as a fxn of Nu and WR mass

#replace int values, fxn name, and abs path name in calculateGenWrScaleFactors_temp.C
eval "sed 's@NULL@@' calculateGenWrScaleFactors_temp.C > macrozero.C"
eval "sed 's@PTHTOTREES@$pathToTrees@g' macrozero.C > macroone.C"
eval "sed 's@QNUM@$q@g' macroone.C > macrotwo.C"
eval "sed 's@MXWRMSS@$maxWrMass@g' macrotwo.C > macrothree.C"
eval "sed 's@INCRMT@$increment@g' macrothree.C > macrofour.C"
eval "sed 's@MNWRMSS@$wrMass@g' macrofour.C > macrofive.C"
eval "sed 's@MNNUMSS@$minNuMass@g' macrofive.C > calculateGenWrScaleFactors.C"

#delete temporary files created with sed
rm macro*.C

####

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
		eval "sed 's@TREEDIR@$pathToTrees@g' run_two.csh > run_three.csh"
		eval "sed 's@GENDIR@$pathToGenEvts@g' run_three.csh > run_four.csh"
		eval "sed 's/MASSNU/$nuMass/g' run_four.csh > runGenAndOfflineAnalysis_WR_${wrMass}_NU_${nuMass}_${q}.csh"
		
		rm run_one.csh run_two.csh run_three.csh run_four.csh

		#replacements in job submission file
		eval "sed 's/NUM/$q/g' genAndRunOfflineAnalysis > job_one"
		eval "sed 's/MMAASS/$wrMass/g' job_one > job_two"
		eval "sed 's@EXEDIR@$pathToExe@g' job_two > job_three"
		eval "sed 's/MASSNU/$nuMass/g' job_three > genAndRunOfflineAnalysis_WR_${wrMass}_NU_${nuMass}_${q}"
		rm job_one job_two job_three

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

