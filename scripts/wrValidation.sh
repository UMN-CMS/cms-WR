#!/bin/bash

#lead and sublead jet pt cuts used to determine dy scaling factors
ljPtCut=40
sjPtCut=40

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

#now make and run analysis.cpp #make lib and outputFromAnalysis directories if they don't already exist
eval "mkdir -p lib"
eval "mkdir -p outputFromAnalysis"
eval "make"

echo "about to run analysis.cpp to produce plotting trees without MLL weights from dytagandprobe minitrees"

tAndPdatasetNames=('data' 'DYMAD' 'DYPOWHEG' 'DYAMC')
for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandPNoMllSf_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandPNoMllSf_${tAndPdatasetNames[$r]}.txt &"
done

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
eval "sed 's@leadJetPtCut = LEADJETPT@leadJetPtCut = $ljPtCut@' <calculateDyScaleFactors_temp.C >tempNoSfOne.C"
eval "sed 's@subLeadJetPtCut = SUBJETPT@subLeadJetPtCut = $sjPtCut@' <tempNoSfOne.C >tempNoSfTwo.C"
eval "sed 's@idealLeadJetPt = 1@idealLeadJetPt = $ljPtCut@' <tempNoSfTwo.C >tempNoSfThree.C"
eval "sed 's@idealSubleadJetPt = 1@idealSubleadJetPt = $sjPtCut@' <tempNoSfThree.C >tempNoSfFour.C"
eval "sed 's@useMllReweighted = true@useMllReweighted = false@' <tempNoSfFour.C >calculateDyScaleFactors.C"
eval "root -l -b runCalculateDyScaleFactors.C"

rm calculateDyScaleFactors.C tempNoSfOne.C tempNoSfTwo.C tempNoSfThree.C tempNoSfFour.C
eval "mkdir -p validationPlots"
eval "mv *.png *.pdf validationPlots/."
eval "git add ../configs/dyScaleFactors.txt"
echo "MLL scale factors have been calculated and written to configs/dyScaleFactors.txt"
eval "cd ../."

echo "about to reprocess the minitrees with dy scaling factors"

#W dataset is WJets
sidebandDatasetNames=('data' 'TT' 'W' 'WZ' 'ZZ' 'DYMAD' 'DYPOWHEG' 'DYAMC')
for r in ${!sidebandDatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors false --isLowDiLepton true --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkMuMuLowDileptonMass_${sidebandDatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors false --isLowDiLepton true --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEELowDileptonMass_${sidebandDatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EMu --ignoreDyScaleFactors false --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEMu_${sidebandDatasetNames[$r]}.txt &"
done

tAndPdatasetNames=('data' 'DYMAD' 'DYPOWHEG' 'DYAMC')
for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandP_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandP_${tAndPdatasetNames[$r]}.txt &"
done

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
eval "./scripts/dyScalingStudy.sh"
eval "cd test/"
eval "mv *.png *.pdf validationPlots/."
eval "cd ../."
eval "git commit -m 'finished running validation and committing ZPeakRegionIntegrals.txt, miniTreeEntries.dat and dyScaleFactors.txt in configs dir' "

