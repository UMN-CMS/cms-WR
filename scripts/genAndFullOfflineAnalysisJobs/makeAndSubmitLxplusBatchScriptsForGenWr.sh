#!/bin/bash

#run this script from the scripts/genAndFullOfflineAnalysis/. directory

#lowest values for WR and Nu masses
#minNuMass and nuMass should be equal
#nuMass is used as an iterator which is frequently changed, while minNuMass is used
#as a reset point for nuMass when wrMass is incremented
#the value of wrMass set here determines the lowest WR mass which will be simulated
wrMass=800
reducedWrMass=600	#should always be 200 less than wrMass
nuMass=25
minNuMass=25
maxWrMass=2400
increment=25
wrIncrement=100

label="genWrToMuMuJJFullOfflineAnalysis_WR"
masterBatchSubDir="batchSubmFilesAndLogDirs"
jobStartingDir="$PWD/../.."
gridProxyPath="/afs/cern.ch/user/s/skalafut/x509up_u38430"
leptonChannel="MuMu"
outputFileDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/analyzedGen/withoutGenNuFilter"	#do not add a fwd slash at the end of this string

eval "mkdir -p $outputFileDir"

#make a directory for the crab python files and crab project directories
eval "mkdir -p $masterBatchSubDir"

#q is a number appended to the end of the output .root file names
#nevts is the number of evts generated at each mass point (default is 15000)
q=1
nevts=10000

#if there are more than maxRunning jobs running, delay the submission of jobs until the number of jobs drops below safeLowerBoundNumRunningJobs
currentSubmittedJobs=0
maxRunning=120

####
##NOTE these 6 sed commands are added to maintain functionality of the root macro
##which calculates the GEN efficiency values as a fxn of Nu and WR mass

#replace int values, fxn name, and abs path name in calculateGenWrScaleFactors_temp.C
eval "sed 's@NULL@@' calculateGenWrScaleFactors_temp.C > macrozero.C"
eval "sed 's@PTHTOTREES@$outputFileDir@g' macrozero.C > macroone.C"
eval "sed 's@QNUM@$q@g' macroone.C > macrotwo.C"
eval "sed 's@MXWRMSS@$maxWrMass@g' macrotwo.C > macrothree.C"
eval "sed 's@INCRMT@$increment@g' macrothree.C > macrofour.C"
eval "sed 's@MNWRMSS@$wrMass@g' macrofour.C > macrofive.C"
eval "sed 's@MNNUMSS@$minNuMass@g' macrofive.C > calculateGenWrScaleFactors.C"

#delete temporary files created with sed
rm macro*.C

####

##main job submission
while [ $wrMass -le $maxWrMass ]; do
	while [ $nuMass -lt $wrMass ]; do
		if [ $nuMass -ge 310 -a $nuMass -le $reducedWrMass ]; then
			let nuMass=nuMass+$increment
			continue
		fi

		#replace MMAASS, MASSNU, EVTCT, and NUM in WR_M-UNDEF_ToLNu_M-UNDEF_GEN.py
		eval "sed 's@MMAASS@$wrMass@g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN.py > genOne.py"
		eval "sed 's@MASSNU@$nuMass@g' genOne.py > genTwo.py"
		eval "sed 's/NUM/$q/g' genTwo.py > genThree.py"
		eval "sed 's@EVTCT@$nevts@' genThree.py > WR_M-${wrMass}_ToLNu_M-${nuMass}_GEN_${q}.py"
		eval "mv WR_M-${wrMass}_ToLNu_M-${nuMass}_GEN_${q}.py ../../."
		rm genOne.py genTwo.py genThree.py 

		#replace CHNL, PROXYPATH, LOCALPATH, OUTPATH, SMPL, NUM, MMAASS, and MASSNU in runBatchJobTemp.sh
		#label has the channel information
		eval "sed 's/NUM/$q/g' runBatchJobTemp.sh > jobOne.sh"
		eval "sed 's/MMAASS/$wrMass/g' jobOne.sh > jobTwo.sh"
		eval "sed 's/MASSNU/$nuMass/g' jobTwo.sh > jobThree.sh"
		eval "sed 's@PROXYPATH@$gridProxyPath@' jobThree.sh > jobFour.sh"
		eval "sed 's@LOCALPATH@$jobStartingDir@' jobFour.sh > jobFive.sh"
		eval "sed 's@CHNL@$leptonChannel@' jobFive.sh > jobSix.sh"
		eval "sed 's@SMPL@$label@' jobSix.sh > jobSeven.sh"
		eval "sed 's@OUTPATH@$outputFileDir@' jobSeven.sh > batchJob_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}.sh"

		rm job*.sh
		mv batchJob_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}.sh $masterBatchSubDir
		#echo "bsub -R 'rusage[mem=1200] && pool>2000' -q 1nh -J analyze_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}_job < batchJob_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}.sh"
		eval "cd $masterBatchSubDir"
		eval "bsub -R 'rusage[mem=1200] && pool>2000' -q 1nh -J analyze_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}_job < batchJob_${q}_${label}_${wrMass}_ToLNu_M_${nuMass}.sh"
		eval "cd ../."
	
		#increment nuMass before restarting the loop
		let nuMass=nuMass+$increment

		#increment currentSubmittedJobs
		let currentSubmittedJobs=currentSubmittedJobs+1

		#if the number of running jobs is large, delay the submission of more jobs
		if [ $currentSubmittedJobs -ge $maxRunning ]; then
			let currentSubmittedJobs=0
			eval "sleep 40m"
		fi

	done
	#reset the nuMass value to the lower bound value before restarting the while loop
	#with a higher wrMass value
	let nuMass=minNuMass
	let wrMass=wrMass+$wrIncrement
	let reducedWrMass=wrMass-200

done

###only for regenerating files and resubmitting jobs which failed
###comment out the lines below when doing main submission defined above
###wrMasses and nuMasses must have the same number of elements
#wrMasses=(1500 1600 2300 2600 2800 3000 3200 3200 3400 3500 3500 3500 3600 3600)
#nuMasses=(100 1400 2000 2300 2600 2900 800 2300 3100 400 2400 3100 400 700)
#
#for i in ${!wrMasses[*]}
#do
#	#replace MMAASS, MASSNU, EVTCT, and NUM in WR_M-UNDEF_ToLNu_M-UNDEF_GEN.py
#	eval "sed 's@MMAASS@${wrMasses[$i]}@g' WR_M-UNDEF_ToLNu_M-UNDEF_GEN.py > genOne.py"
#	eval "sed 's@MASSNU@${nuMasses[$i]}@g' genOne.py > genTwo.py"
#	eval "sed 's/NUM/$q/g' genTwo.py > genThree.py"
#	eval "sed 's@EVTCT@$nevts@' genThree.py > WR_M-${wrMasses[$i]}_ToLNu_M-${nuMasses[$i]}_GEN_${q}.py"
#	eval "mv WR_M-${wrMasses[$i]}_ToLNu_M-${nuMasses[$i]}_GEN_${q}.py ../../."
#	rm genOne.py genTwo.py genThree.py 
#
#	#replace CHNL, PROXYPATH, LOCALPATH, OUTPATH, SMPL, NUM, MMAASS, and MASSNU in runBatchJobTemp.sh
#	#label has the channel information
#	eval "sed 's/NUM/$q/g' runBatchJobTemp.sh > jobOne.sh"
#	eval "sed 's/MMAASS/${wrMasses[$i]}/g' jobOne.sh > jobTwo.sh"
#	eval "sed 's/MASSNU/${nuMasses[$i]}/g' jobTwo.sh > jobThree.sh"
#	eval "sed 's@PROXYPATH@$gridProxyPath@' jobThree.sh > jobFour.sh"
#	eval "sed 's@LOCALPATH@$jobStartingDir@' jobFour.sh > jobFive.sh"
#	eval "sed 's@CHNL@$leptonChannel@' jobFive.sh > jobSix.sh"
#	eval "sed 's@SMPL@$label@' jobSix.sh > jobSeven.sh"
#	eval "sed 's@OUTPATH@$outputFileDir@' jobSeven.sh > batchJob_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}.sh"
#
#	rm job*.sh
#	mv batchJob_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}.sh $masterBatchSubDir
#	eval "cd $masterBatchSubDir"
#	#echo "bsub -R 'rusage[mem=1200] && pool>2000' -q 1nh -J analyze_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}_job < batchJob_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}.sh"
#	eval "bsub -R 'rusage[mem=1200] && pool>2000' -q 1nh -J analyze_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}_job < batchJob_${q}_${label}_${wrMasses[$i]}_ToLNu_M_${nuMasses[$i]}.sh"
#	eval "cd ../."
#
#done
