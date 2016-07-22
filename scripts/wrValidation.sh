#!/bin/bash
###this script can be run at any time from the cmsWR/. directory as long as minitrees are available
###this script processes minitrees from collision data and MC datasets with analysis.cpp, and uses the resulting root files
###to make plots which compare data to MC in several sidebands and channels
###the only parameters which should be modified in this script are the four cut values ljPtCut, sjPtCut, llPtCut, and slPtCut
###defined below. ljPtCut (leading jet) and sjPtCut (subleading jet) define the leading and subleading jet pt cuts used in
###dytagandprobe sideband events, and llPtcut (leading lepton) and slPtCut (subleading lepton) define the leading and
###subleading lepton pt cuts used in dytagandprobe sideband events. Specifically, these are the cut values for which
###dy MC scale factors will be saved to dyScaleFactors.txt. The lepton pt cuts must agree with the lepton pt cuts used in
###the isPassingLooseCuts() method defined in the Selector class.


#lead and sublead jet pt cuts used to determine dy scaling factors
ljPtCut=40
sjPtCut=40
llPtCut=35
slPtCut=35

tAndPdatasetNames=('data' 'DYPOWHEG' 'DYAMC' 'DYMADHT')
#W is WJetsToLNu
sidebandDatasetNames=('data' 'TT' 'W' 'WZ' 'ZZ' 'DYMADHT' 'DYPOWHEG' 'DYAMC')


#cd to cmsWR/. and run the python script test/miniTreeDumpCheck.py to check the number of entries in minitrees
eval "rm -r test/validationPlots/"

eval "echo -n '#' > data/2015-v1/miniTreeEntries.dat"
eval "echo priorCommitNumberForWhichTheseEventCountNumbersAreValid >> data/2015-v1/miniTreeEntries.dat"
eval "echo -n '#' >> data/2015-v1/miniTreeEntries.dat"
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> data/2015-v1/miniTreeEntries.dat"

#this next line removes the bizarre ^[[?1034h code which is written to the first line of a file
#when bash is used to redirect the output of a python script to a file
eval "TERM=linux"
eval 'python test/miniTreeDumpCheck.py >> data/2015-v1/miniTreeEntries.dat'

echo "finished checking num entries in minitrees"

#now make and run analysis.cpp
#make lib and outputFromAnalysis directories if they don't already exist
eval "mkdir -p lib"
eval "mkdir -p outputFromAnalysis"
eval "make"

echo "about to run analysis.cpp to produce plotting trees without MLL scale factor from dytagandprobe minitrees"

#run this to allow DY MLL scale factors to be calculated
for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandPWithoutMllSF_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors true --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandPWithoutMllSF_${tAndPdatasetNames[$r]}.txt &"
done

wait
#delay this script until all of the dytagandprobe analysis processes are finished, then calculate the dy scaling factors
#while :
#do
#	x=1
#	if pgrep "[a]nalysis" > /dev/null
#	then
#		y=1
#	else
#		break
#	fi
#done
echo "all dytagandprobe analysis processes are finished"

echo "about to process TTBar MC with double ele and double mu signal region requirements"
#process the TTBar dataset with EEJJ and MuMuJJ signal region requirements to then calculate ttbar scale factors (ee to emu and mumu to emu)
eval "./bin/analysis -c MuMu -m TT >& outputFromAnalysis/checkMuMuSignalWithSFs_TT.txt &"
eval "./bin/analysis -c EE -m TT >& outputFromAnalysis/checkEESignalWithSFs_TT.txt &"

wait

echo "all ttbar signal region analysis processes are finished"
echo "starting to calculate dy scaling factors"
eval "cd test"

#the jet pt cut applied to both jets in calculateDyScaleFactors.C
jetCutVals=('-10.0' '10.0' '20.0' '30.0' '40.0' '50.0')

#the lepton pt cut values applied to the leading and subleading lepton in calculateDyScaleFactors.C
#leadLeptCutVals and subleadLeptCutVals must have the same number of elements
#events are selected by analysis.cpp that have two leptons with pt>=35, so pt cuts below 35 are pointless
leadLeptCutVals=('40.0' '40.0')
subleadLeptCutVals=('35.0' '40.0')

#overwrite data/2015-v1/allDyScaleFactors.txt if it already exists
eval "echo 'channel   DY   data/DY   leadLeptPtCut   subleadLeptPtCut   leadJetPtCut   subleadJetPtCut' > ../data/2015-v1/allDyScaleFactors.txt"
for e in ${!jetCutVals[*]}
do
	#leadJetPtCut and subLeadJetPtCut are the actual jet pt cuts which are applied when histograms are filled
	#these values need to be varied
	eval "sed 's@leadJetPtCut = LEADJETPT@leadJetPtCut = ${jetCutVals[$e]}@' <calculateDyScaleFactors_temp.C >tempNoSfOne.C"
	eval "sed 's@subLeadJetPtCut = SUBJETPT@subLeadJetPtCut = ${jetCutVals[$e]}@' <tempNoSfOne.C >tempNoSfTwo.C"
	
	#idealLeadJetPt and idealSubleadJetPt determine which txt file in the data/2015-v1 dir the data over DY ratios are written
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

#make directory where all plots and text files will be stored
eval "mkdir -p validationPlots"
eval "mv *.png *.pdf validationPlots/."
eval "git add ../data/2015-v1/dyScaleFactors.txt ../data/2015-v1/allDyScaleFactors.txt"
echo "MLL scale factors have been calculated and written to data/2015-v1/dyScaleFactors.txt and data/2015-v1/allDyScaleFactors.txt"
eval "cd ../."

echo "about to reprocess the minitrees with dy scaling factors"

#W dataset is WJets
for r in ${!sidebandDatasetNames[*]}
do
	eval "./bin/analysis -c MuMu --ignoreDyScaleFactors false --isLowDiLepton true -m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkMuMuLowDileptonMass_${sidebandDatasetNames[$r]}_withMllSf.txt &"
	eval "./bin/analysis -c EE --ignoreDyScaleFactors false --isLowDiLepton true -m ${sidebandDatasetNames[$r]} >& outputFromAnalysis/checkEELowDileptonMass_${sidebandDatasetNames[$r]}_withMllSf.txt &" 
	eval "./bin/analysis -c EMu --ignoreDyScaleFactors false -m ${sidebandDatasetNames[$r]} --cut_channel EE >& outputFromAnalysis/checkEMu_${sidebandDatasetNames[$r]}_withMllSf_EEcutChannel.txt &"
done

wait
#wait until all low dilepton mass and flavor sideband processing jobs are finished before starting dytagandprobe processing. All events in each minitree will be loaded into
#memory, and this requires a large amount of memory for the dytagandprobe minitrees, the dyamc minitrees in particular
echo "all low dilepton mass and flavor sideband analysis processes using MLL scale factors are finished"

for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandPWithMllSF_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --ignoreDyScaleFactors false --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandPWithMllSF_${tAndPdatasetNames[$r]}.txt &"
done

wait 
echo "all dytagandprobe analysis processes using MLL scale factors are finished"

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

eval "cp ../data/2015-v1/miniTreeEntries.dat validationPlots/."
eval "cp ../data/2015-v1/dyScaleFactors.txt validationPlots/."
eval "root -l -b runCombinedMiniPlotterForDYTandPMuMu.C"
eval "root -l -b runCombinedMiniPlotterForDYTandPEE.C"
eval "root -l -b runCombinedMiniPlotterMuMu.C"
eval "root -l -b runCombinedMiniPlotterEE.C"
eval "root -l -b runCombinedMiniPlotterEMu.C"
eval "cp validationPlots/*.txt ../data/2015-v1/."
eval "git add ../data/2015-v1/ZPeakRegionIntegrals.txt ../data/2015-v1/miniTreeEntries.dat"
rm runCombinedMiniPlotterForDYTandPEE.C runCombinedMiniPlotterEE.C runCombinedMiniPlotterEMu.C
rm combinedMiniPlotterForDYTandPEE.C combinedMiniPlotterEE.C combinedMiniPlotterEMu.C tempEMu.C
rm combinedMiniPlotterForDYTandPEE_C*
rm combinedMiniPlotterEE_C*
rm combinedMiniPlotterEMu_C*


eval "cd ../."
#dyScalingStudy.sh should be run after the scaling factors are calculated and applied
#running this script is time consuming, and is not needed as part of the minitree validation framework
#eval "./scripts/dyScalingStudy.sh"
#eval "cd test/"
#eval "mv *.png *.pdf validationPlots/."
#eval "cd ../."

#commit the text files which summarize the status of the minitrees and data to DY scale factors
eval "git commit -m 'finished running validation and committing ZPeakRegionIntegrals.txt, miniTreeEntries.dat and dyScaleFactors.txt in data/2015-v1 dir' "

