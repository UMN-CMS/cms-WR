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
label="genWrToEEJJFullOfflineAnalysis_WR"
masterCrabDir="crabPythonFilesAndProjDirs"

#make a directory for the crab python files and crab project directories
eval "mkdir -p $masterCrabDir"

#q is a number appended to the end of the output .root file names
q=1

#absolute path names to output file directories
#where GEN events and TTree files will be stored
pathToTrees='/uscms/home/skalafut/nobackup/WR_starting2015/genWrTTrees'

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

		#replace SMPL, NUM, MMAASS, and MASSNU in analyze_privateWrGen_crab.py
		#label has the channel information, so ToLNu is fine
		eval "sed 's/NUM/$q/g' analyze_privateWrGen_crab.py > crabOne.py"
		eval "sed 's/MMAASS/$wrMass/g' crabOne.py > crabTwo.py"
		eval "sed 's/MASSNU/$nuMass/g' crabTwo.py > crabThree.py"
		eval "sed 's/SMPL/$label/g' crabThree.py > analyze_${label}_${wrMass}_ToLNu_M-${nuMass}_crab_${q}.py"
		rm crabOne.py crabTwo.py crabThree.py
		mv analyze_${label}_${wrMass}_ToLNu_M-${nuMass}_crab_${q}.py $masterCrabDir
		eval "cd $masterCrabDir"
		eval "crab submit -c analyze_${label}_${wrMass}_ToLNu_M-${nuMass}_crab_${q}.py"
		eval "cd ../."
	
		#increment nuMass before restarting the loop
		let nuMass=nuMass+$increment
	done
	#reset the nuMass value to the lower bound value before restarting the while loop
	#with a higher wrMass value
	let nuMass=minNuMass
	let wrMass=wrMass+$increment

done

