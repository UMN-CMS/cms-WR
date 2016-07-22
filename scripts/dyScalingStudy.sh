#!/bin/bash

###only use this script AFTER the dy scale factors have been calculated and used by analysis.cpp to reprocess the dytagandprobe minitrees for DY MC and data
###this script studies the dytagandprobe data and MC after the DY mll scale factors have been applied
###as a function of jet pt cuts
###this script can be run by itself, or within wrValidation.sh
###the only parameters which can be changed without consulting the lepton and jet default pts and pt cuts in the Selector class and isPassingLooseCuts() are
###maxJetPt and increment. initialPtCut, minLeadJetPt, and minSubleadJetPt must be equal.

leadLeptPtCut=35
subleadLeptPtCut=35
initialPtCut=-10
maxJetPt=40
minLeadJetPt=-10
increment=10
minSubleadJetPt=-10
lowestPositiveCutVal=0
switchPoint=-1

eval "cd test"
#make sure that the Mll reweighted trees are used with dR jet lepton separation required
#and the lead and sublead lepton pt cuts are set to appropriate values
eval "sed 's@useMllReweighted = false@useMllReweighted = true@' <calculateDyScaleFactors_temp.C >tempB.C"
eval "sed 's@LEADLEPTPT@$leadLeptPtCut@' <tempB.C >tempC.C"
eval "sed 's@SUBLEPTPT@$subleadLeptPtCut@' <tempC.C >tempA.C"
eval "sed 's@requireSeparatedLeptonsAndJets = false@requireSeparatedLeptonsAndJets = true@' <tempA.C >tempOne.C"
rm tempA.C tempB.C tempC.C

#make one set of plots with the requirement of two jets with pt above 0, and no jet lepton dR separation
#then proceed to make many sets of plots with various jet pt requirements with jet lepton dR separation
eval "sed 's@requireSeparatedLeptonsAndJets = true@requireSeparatedLeptonsAndJets = false@' <tempOne.C >tempNoDrOne.C"
eval "sed 's@LEADJETPT@$lowestPositiveCutVal@' <tempNoDrOne.C >tempNoDrTwo.C"
eval "sed 's@SUBJETPT@$lowestPositiveCutVal@' <tempNoDrTwo.C >calculateDyScaleFactors.C"
eval "root -l -b runCalculateDyScaleFactors.C"
rm tempNoDrOne.C tempNoDrTwo.C calculateDyScaleFactors.C

#now make many sets of plots with various jet pt requirements with jet lepton dR separation required
while [ $minLeadJetPt -le $maxJetPt ]; do
	while [ $minSubleadJetPt -le $minLeadJetPt ]; do
		#set the jet pt cuts to the current values
		eval "sed 's@LEADJETPT@$minLeadJetPt@' <tempOne.C >tempTwo.C"
		eval "sed 's@SUBJETPT@$minSubleadJetPt@' <tempTwo.C >calculateDyScaleFactors.C"
		
		#run the macro calculateDyScaleFactors.C and produce the plots
		eval "root -l -b runCalculateDyScaleFactors.C"

		#remove the files which will be overwritten by the next loop iteration
		rm tempTwo.C calculateDyScaleFactors.C
	
		#increase the sublead jet pt cut
		if [ $minSubleadJetPt -gt $switchPoint ]; then
			#increase the sublead jet pt cut by increment if the cut is already greater than zero
			let minSubleadJetPt=($minSubleadJetPt+$increment)
		fi
		
		if [ $minSubleadJetPt -lt $switchPoint ]; then
			#if minSubleadJetPt is negative, then increase it to a nonnegative value
			let minSubleadJetPt=$lowestPositiveCutVal
		fi
	done
	#reset the sublead jet pt requirement and increase the lead jet pt cut
	let minSubleadJetPt=$initialPtCut
	if [ $minLeadJetPt -gt $switchPoint ]; then
		#increase the lead jet pt cut by increment if the cut is already greater than zero
		let minLeadJetPt=($minLeadJetPt+$increment)
	fi

	if [ $minLeadJetPt -lt $switchPoint ]; then
		#if minLeadJetPt is negative, then increase it to a nonnegative value
		let minLeadJetPt=$lowestPositiveCutVal
	fi
done

rm tempOne.C

