#!/bin/bash

#bkgndProcess=ttBar
bkgndProcess=dyPlusJets

#q btwn 0 and 275 for ttbar
#q btwn 289 and 501 for dyPlusJets
for q in {289..501} 
do
	#request_disk in KB, request_memory in MB
	#echo 'condor_submit request_disk=8000000 request_memory=5000 recoAnalysis_${bkgndProcess}_$q'
	eval 'condor_submit request_disk=8000000 request_memory=5000 recoAnalysis_${bkgndProcess}_$q'

done


