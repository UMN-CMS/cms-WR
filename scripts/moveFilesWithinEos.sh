#!/bin/bash
#execute this in the ExoAnalysis/cmsWR/ directory

#source this file to read the correct datasets .dat file
source configs/2015-v1.conf

eosSkimPath='/eos/cms/store/group/phys_exotica/leptonsPlusJets/WR/skims/'
#the awk in the next line is not needed, but it is added just for reference. this next line is based on simplerMakeAndRunCrabSkimScripts.sh
topLevelDirs=(`eos ls $eosSkimPath | grep pythia | awk '{print $1}'`)

#datasetFile defined in configs dir, 2015-v1.conf
newDirsToMake=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)

##only need to run this once
##first make the new directories will the output files will be copied
#for j in ${!newDirsToMake[*]}
#do
#	#make a new directory with the shortened dataset name plus skimProductionTAG defined in 2015-v1.conf
#	#eos mkdir -p $eosSkimPath${newDirsToMake[$j]}$skimProductionTAG
#	echo "eos mkdir -p $eosSkimPath${newDirsToMake[$j]}$skimProductionTAG"
#
#done

#cpTag='root://eoscms.cern.ch/'
#for i in ${!topLevelDirs[*]}
#do
#	#get the list of subdirectories and files within one skim directory, like WW or QCD EMEnriched pt30to50
#	#element i in topLevelDirs corresponds to element i in newDirsToMake, which have already been made
#	#all elements in skimFiles are root files which were selected by the skims
#	skimFiles=(`eos find $eosSkimPath${topLevelDirs[$i]} | grep output`)
#
#	for k in ${!skimFiles[*]}
#	do
#		let y=k+1
#		#RUN this echo command and comment out code further below before trying to actually cp any files between dirs on eos.
#		#make sure the dir name in the skimFiles element matches with the element in newDirsToMake, otherwise
#		#a skim from one dataset will end up in a skim of a different dataset
#		#echo "xrdcp $cpTag${skimFiles[$k]} ."
#		#echo "xrdcp output_${y}.root $cpTag$eosSkimPath${newDirsToMake[$i]}$skimProductionTAG/"
#		#echo "rm output_${y}.root"
#
#		#copy each element from the original dir to the local dir, and the local dir to the new dir
#		eval "xrdcp $cpTag${skimFiles[$k]} ."
#		eval "xrdcp output_${y}.root $cpTag$eosSkimPath${newDirsToMake[$i]}$skimProductionTAG/"
#		eval "rm output_${y}.root"
#
#	done
#
#done

#remove the original directory on eos where the skim output files were stored
#for q in ${!topLevelDirs[*]}
#do
#	#RUN the echo command before actually removing any directory
#	#echo "eos rm -r $eosSkimPath${topLevelDirs[$q]}"
#	
#	#eos rm -r $eosSkimPath${topLevelDirs[$q]}
#
#done


##only need to run this once
##make new directories in WR/tuples/ where the minitrees will be stored
eosTuplePath='/eos/cms/store/group/phys_exotica/leptonsPlusJets/WR/tuples/'
for j in ${!newDirsToMake[*]}
do
	#make a new directory with the shortened dataset name plus productionTAG defined in 2015-v1.conf
	eos mkdir -p $eosTuplePath${newDirsToMake[$j]}$productionTAG
	#echo "eos mkdir -p $eosTuplePath${newDirsToMake[$j]}$productionTAG"

done




