#!/bin/bash

##cd to cmsWR/. and run the python script test/miniTreeDumpCheck.py to check the number of entries in minitrees
eval "rm -r test/validationPlots/"

echo -n '#' > configs/miniTreeEntries.dat
echo priorCommitNumberForWhichTheseEventCountNumbersAreValid >> configs/miniTreeEntries.dat
echo -n '#' >> configs/miniTreeEntries.dat
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> configs/miniTreeEntries.dat"

#this next line removes the bizarre ^[[?1034h code which is written to the first line of a file
#when bash is used to redirect the output of a python script to a file
eval "TERM=linux"
eval 'python test/miniTreeDumpCheck.py >> configs/miniTreeEntries.dat'

#now make and run analysis.cpp #make lib and outputFromAnalysis directories if they don't already exist
eval "mkdir -p lib"
eval "mkdir -p outputFromAnalysis"
eval "make"

#W dataset is WJets
sidebandDatasetNames=('data' 'TT' 'W' 'WZ' 'ZZ' 'DYMAD')
for r in ${!sidebandDatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --isLowDiLepton true --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkMuMuLowDileptonMass_${sidebandDatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --isLowDiLepton true --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEELowDileptonMass_${sidebandDatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EMu --m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEMu_${sidebandDatasetNames[$r]}.txt &"
done

tAndPdatasetNames=('data' 'DYMAD' 'DYPOWHEG' 'DYAMC')
for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandP_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandP_${tAndPdatasetNames[$r]}.txt &"
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
eval "mkdir -p validationPlots"

#overwrite ZPeakRegionIntegrals.txt if it already exists
echo -n '#' > validationPlots/ZPeakRegionIntegrals.txt
echo priorCommitNumberForWhichTheseEventCountNumbersAreValid >> validationPlots/ZPeakRegionIntegrals.txt
echo -n '#' >> validationPlots/ZPeakRegionIntegrals.txt
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> validationPlots/ZPeakRegionIntegrals.txt"

eval "cp ../configs/miniTreeEntries.dat validationPlots/."
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





