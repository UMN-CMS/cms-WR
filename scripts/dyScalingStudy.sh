#!/bin/bash

#NOT WORKING YET DO NOT USE
#NOT WORKING YET DO NOT USE

#only use this script AFTER the dy scale factors have been calculated and used to reprocess the dytagandprobe minitrees for DY MC and data with analysis.cpp
initialPtCut=-10
maxJetPt=40
minLeadJetPt=-10
increment=10
minSubleadJetPt=-10
lowestPositiveCutVal=10

while [$minLeadJetPt -le $maxJetPt]; do
	while [$minSubleadJetPt -le $minLeadJetPt]; do
		#make plots with the existing pt cuts
		cd test
		eval "sed 's@@@' > "

		#now increase the sublead jet pt cut
		if [$minSubleadJetPt > $initialPtCut]; then
			#increase the sublead jet pt cut by increment if the cut is already greater than initialPtCut
			let minSubleadJetPt=minSubleadJetPt+$increment
		fi
		
		if [$minSubleadJetPt == $initialPtCut]; then
			#if minSubleadJetPt is negative, then increase it to a nonnegative value
			let minSubleadJetPt=$lowestPositiveCutVal
		fi
		#let minSubleadJetPt=minSubleadJetPt+$increment
	done
	#reset the sublead jet pt requirement and increase the lead jet pt cut by increment
	let minSubleadJetPt=initialPtCut
	let minLeadJetPt=minLeadJetPt+$increment
done

