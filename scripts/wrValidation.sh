#!/bin/bash

#lead and sublead jet pt cuts used to determine dy scaling factors
ljPtCut=40
sjPtCut=40
llPtCut=35
slPtCut=35

tAndPdatasetNames=('data' 'DYPOWHEG' 'DYAMC' 'DYMADHT')
sidebandDatasetNames=('data' 'TT' 'W' 'WZ' 'ZZ' 'DYMADHT' 'DYPOWHEG' 'DYAMC')

#cd to cmsWR/. and run the python script test/miniTreeDumpCheck.py to check the number of entries in minitrees
eval "rm -r test/validationPlots/"

echo -n '#' > configs/miniTreeEntries.dat
echo priorCommitNumberForWhichTheseEventCountNumbersAreValid >> configs/miniTreeEntries.dat
echo -n '#' >> configs/miniTreeEntries.dat
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> configs/miniTreeEntries.dat"

#this next line removes the bizarre ^[[?1034h code which is written to the first line of a file
#when bash is used to redirect the output of a python script to a file
eval "TERM=linux"
eval 'python test/miniTreeDumpCheck.py >> configs/miniTreeEntries.dat'

echo "finished checking num entries in minitrees"

now make and run analysis.cpp #make lib and outputFromAnalysis directories if they don't already exist
eval "mkdir -p lib"
eval "mkdir -p outputFromAnalysis"
eval "make"

echo "about to run analysis.cpp to produce plotting trees without MLL weights from dytagandprobe minitrees"

for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandPNoMllSf_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandPNoMllSf_${tAndPdatasetNames[$r]}.txt &"
done

wait
#delay this script until all of the dytagandprobe analysis processes are finished, then calculate the dy scaling factors
while :
do
	x=1
	if pgrep "[a]nalysis" > /dev/null
	then
		y=1
	else
		break
	fi
done

echo "all dytagandprobe analysis processes are finished   starting to calculate dy scaling factors"
eval "cd test"

#this code replaces the functionality of the dyScalingStudy.sh before MC is rescaled to data
#the jet pt cut applied to both jets in calculateDyScaleFactors.C
jetCutVals=('-10.0' '10.0' '30.0' '35.0' '40.0' '45.0')

#the lepton pt cut values applied to the leading and subleading lepton in calculateDyScaleFactors.C
#leadLeptCutVals and subleadLeptCutVals must have the same number of elements
leadLeptCutVals=('33.0' '33.0' '33.0' '35.0' '35.0' '35.0' '37.0' '37.0' '40.0' '40.0')
subleadLeptCutVals=('20.0' '25.0' '29.0' '25.0' '29.0' '31.0' '29.0' '33.0' '33.0' '35.0')

#overwrite configs/allDyScaleFactors.txt if it already exists
echo 'channel   DY   data/DY   leadLeptPtCut   subleadLeptPtCut   leadJetPtCut   subleadJetPtCut' > ../configs/allDyScaleFactors.txt
for e in ${!jetCutVals[*]}
do
	#leadJetPtCut and subLeadJetPtCut are the actual jet pt cuts which are applied when histograms are filled
	#these values need to be varied
	eval "sed 's@leadJetPtCut = LEADJETPT@leadJetPtCut = ${jetCutVals[$e]}@' <calculateDyScaleFactors_temp.C >tempNoSfOne.C"
	eval "sed 's@subLeadJetPtCut = SUBJETPT@subLeadJetPtCut = ${jetCutVals[$e]}@' <tempNoSfOne.C >tempNoSfTwo.C"
	
	#idealLeadJetPt and idealSubleadJetPt determine which txt file in the configs dir the data over DY ratios are written
	#these values do not need to be varied
	eval "sed 's@idealLeadJetPt = 1@idealLeadJetPt = $ljPtCut@' <tempNoSfTwo.C >tempNoSfThree.C"
	eval "sed 's@idealSubleadJetPt = 1@idealSubleadJetPt = $sjPtCut@' <tempNoSfThree.C >tempNoSfFour.C"
	eval "sed 's@useMllReweighted = true@useMllReweighted = false@' <tempNoSfFour.C >tempNoSfFive.C"

	for l in ${!leadLeptCutVals[*]}
	do
		eval "sed 's@LEADLEPTPT@${leadLeptCutVals[$l]}@' <tempNoSfFive.C >tempNoSfSix.C"
		eval "sed 's@SUBLEPTPT@${subleadLeptCutVals[$l]}@' <tempNoSfSix.C >calculateDyScaleFactors.C"
		eval "root -l -b runCalculateDyScaleFactors.C"
		rm calculateDyScaleFactors.C tempNoSfSix.C
	done
	
	rm tempNoSf*.C
	
done

eval "mkdir -p validationPlots"
eval "mv *.png *.pdf validationPlots/."
eval "git add ../configs/dyScaleFactors.txt ../configs/allDyScaleFactors.txt"
echo "MLL scale factors have been calculated and written to configs/dyScaleFactors.txt and configs/allDyScaleFactors.txt"
eval "cd ../."

echo "about to reprocess the minitrees with dy scaling factors"

##W dataset is WJets
for r in ${!sidebandDatasetNames[*]}
do
	eval "./bin/analysis -c MuMu --ignoreDyScaleFactors false --isLowDiLepton true -m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkMuMuLowDileptonMass_${sidebandDatasetNames[$r]}.txt &"
	eval "./bin/analysis -c EE --ignoreDyScaleFactors false --isLowDiLepton true -m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEELowDileptonMass_${sidebandDatasetNames[$r]}.txt &" 
	eval "./bin/analysis -c EMu --ignoreDyScaleFactors false -m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEMu_${sidebandDatasetNames[$r]}.txt &"
done

for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandP_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandP_${tAndPdatasetNames[$r]}.txt &"
done

wait 
#delay this script until all of the analysis processes are finished
while :
do
	x=1
	if pgrep "[a]nalysis" > /dev/null
	then
		y=1
	else
		break
	fi
done

echo "all analysis processes are finished"

#now make duplicates of the plotting macros to make MuMu, EE, and EMu plots
#and run the plotting macros
#one plotting macro is used for DYTagAndProbe, which makes plots with data overlaid onto multiple curves from different DY MC samples
#the other plotting macro is used for the flavor and low dilepton mass sidebands, which show data overlaid on STACKED backgrounds
eval "cd test"
eval "sed 's@MuMu@EE@g' combinedMiniPlotterForDYTandPMuMu.C > combinedMiniPlotterForDYTandPEE.C"
eval "sed 's@MuMu@EE@g' runCombinedMiniPlotterForDYTandPMuMu.C > runCombinedMiniPlotterForDYTandPEE.C"
eval "sed 's@MuMu@EE@g' combinedMiniPlotterMuMu.C > combinedMiniPlotterEE.C"
eval "sed 's@MuMu@EE@g' runCombinedMiniPlotterMuMu.C > runCombinedMiniPlotterEE.C"
eval "sed 's@MuMu@EMu@g' combinedMiniPlotterMuMu.C > tempEMu.C"
eval "sed 's@lowdilepton@flavour@g' tempEMu.C > combinedMiniPlotterEMu.C"
eval "sed 's@MuMu@EMu@g' runCombinedMiniPlotterMuMu.C > runCombinedMiniPlotterEMu.C"

#overwrite ZPeakRegionIntegrals.txt if it already exists
echo -n '#' > validationPlots/ZPeakRegionIntegrals.txt
echo priorCommitNumberForWhichTheseEventCountNumbersAreValid >> validationPlots/ZPeakRegionIntegrals.txt
echo -n '#' >> validationPlots/ZPeakRegionIntegrals.txt
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> validationPlots/ZPeakRegionIntegrals.txt"

eval "cp ../configs/miniTreeEntries.dat validationPlots/."
eval "cp ../configs/dyScaleFactors.txt validationPlots/."
eval "root -l -b runCombinedMiniPlotterForDYTandPMuMu.C"
eval "root -l -b runCombinedMiniPlotterForDYTandPEE.C"
eval "root -l -b runCombinedMiniPlotterMuMu.C"
eval "root -l -b runCombinedMiniPlotterEE.C"
eval "root -l -b runCombinedMiniPlotterEMu.C"
eval "cp validationPlots/*.txt ../configs/."
eval "git add ../configs/ZPeakRegionIntegrals.txt ../configs/miniTreeEntries.dat"
rm runCombinedMiniPlotterForDYTandPEE.C runCombinedMiniPlotterEE.C runCombinedMiniPlotterEMu.C
rm combinedMiniPlotterForDYTandPEE.C combinedMiniPlotterEE.C combinedMiniPlotterEMu.C tempEMu.C
rm combinedMiniPlotterForDYTandPEE_C*
rm combinedMiniPlotterEE_C*
rm combinedMiniPlotterEMu_C*

eval "cd ../."
#dyScalingStudy.sh should be run after the scaling factors are calculated and applied
eval "./scripts/dyScalingStudy.sh"
eval "cd test/"
eval "mv *.png *.pdf validationPlots/."
eval "cd ../."
eval "git commit -m 'finished running validation and committing ZPeakRegionIntegrals.txt, miniTreeEntries.dat and dyScaleFactors.txt in configs dir' "

