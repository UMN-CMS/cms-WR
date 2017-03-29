#!/bin/bash
#run this from cmsWR/.

#import productionTAG and datasetFile
source configs/2015-v1.conf

#all input files to process
##for all datasets in configs/missing_datasets.dat
##first get all the dataset short names, and save them to a list

#for bkgnd minitrees
#mcIdentifier=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)

#for WR minitrees
mcIdentifier=(`cat $datasetFile | grep -v '#' | grep 'WR' | awk '{print $1}'`)

eosTuplePath='/eos/cms/store/group/phys_exotica/leptonsPlusJets/WR/tuples/'
eosReadingTag='root://eoscms.cern.ch/'

#now loop over all datasets
for j in ${!mcIdentifier[*]}
do
	#eos mkdir -p $eosTuplePath${mcIdentifier[$j]}$productionTAG
	#only need hadd if there are multiple skim files from one dataset
	#hadd -k -O unmerged-allRange.root ${mcIdentifier[$j]}Part*.root
	cp ${mcIdentifier[$j]}Part1.root unmerged-allRange.root
	xrdcp unmerged-allRange.root $eosReadingTag$eosTuplePath${mcIdentifier[$j]}$productionTAG/.
	rm unmerged-allRange.root
	#echo "hadd -k -O unmerged-allRange.root ${mcIdentifier[$j]}Part*.root"
	#echo "xrdcp unmerged-allRange.root $eosReadingTag$eosTuplePath${mcIdentifier[$j]}$productionTAG/."
	#echo "rm unmerged-allRange.root"

done

