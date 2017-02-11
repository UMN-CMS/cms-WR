#!/bin/bash
#RUN THIS in cmsWR/ dir, not scripts dir

#limit_window_scan.py and similar clones made by makeAndSubmitMassWindowScan.sh calculate the
#expected limit for a given WR mass hypothesis, considering DY and top (top+W and TTBar)
#backgrounds, and print the following string to stdout or a txt file when the limit is calculated
#with a new mass window:
#
# WRmass Numass windowLowerBound windowUpperBound sigEvts bkgndOneEvts bkgndTwoEvts limMinusTwoSigma limMinusOneSigma medianLim limPlusOneSigma limPlusTwoSigma
#
#this script parses txt files which contain lines of strings, in the format of the line above, and
#finds the mass window (lower and upper bound) which minimizes the expected cross sxn limit

chnl='ee'
CHNL='EE'

#make a file which will hold the new mass windows
eval "echo -n '#' > newMassWindows.txt"
eval "echo chnl MWR low hi  >> newMassWindows.txt"

#get the list of files which contain the mass windows and limits
inputFiles=(`ls|grep limitResul|grep -v '1000'|grep -v '1200'|grep $chnl`)

for i in ${!inputFiles[*]}
do
	#for each limit result file, make four lists. one containing the WR mass (constant for every entry),
	#one for the window lower bound, one for the window upper bound, and one for the expected limit
	massHypothesis=(`cat ${inputFiles[$i]}|grep -v 'Profile' |awk '{print $1}'`)
	lowerBound=(`cat ${inputFiles[$i]}|grep -v 'Profile' |awk '{print $3}'`)
	upperBound=(`cat ${inputFiles[$i]}|grep -v 'Profile' |awk '{print $4}'`)
	limit=(`cat ${inputFiles[$i]}|grep -v 'Profile' |awk '{print $10}'`)
	
	#now go through the limit list and find the lowest value
	bestLimitIndex=-1
	bestLimit=185.1
	for n in ${!limit[*]}
	do
		currLim=${limit[$n]}
		compResult=`echo $bestLimit '>' $currLim | bc -l`
		#compResult is 1 if bestLimit is higher than currLim
		if [ $compResult -eq 1 ]; then
			bestLimit=${limit[$n]}
			bestLimitIndex=$n
		fi

	done
	#end loop over list of limits at a given WR mass hypothesis

	#now that the best limit has been found, use bestLimitIndex to print the WR mass, channel, the window lower and upper bound
	#MuMu MWR low hi
	eval "echo $CHNL ${massHypothesis[$bestLimitIndex]} ${lowerBound[$bestLimitIndex]} ${upperBound[$bestLimitIndex]} >> newMassWindows.txt"

done


