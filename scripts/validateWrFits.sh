#!/bin/bash

#run this script starting in in cmsWR/.
eval "cd test/"
echo "now in cmsWR/test/ dir"

#modify Jorge's combinedMiniFitPlot.C to make RooFit plots of RooDataSets and exponential fits applied to those datasets
#EMu fits
eval "sed 's@RGN@flavoursidebandEMu@g' combinedMiniFitPlot_temp.C > tempEMu.C"

eval "sed 's@MODE@data@g' tempEMu.C > combinedMiniFitPlot.C"
eval "root -l -b -q combinedMiniFitPlot.C"
rm combinedMiniFitPlot.C

eval "sed 's@MODE@TT@g' tempEMu.C > combinedMiniFitPlot.C"
eval "root -l -b -q combinedMiniFitPlot.C"
rm combinedMiniFitPlot.C
rm temp*.C 

#EE and MuMu fits
eval "sed 's@RGN@signal_eeEE_withMllWeight@g' combinedMiniFitPlot_temp.C > tempEE.C"
eval "sed 's@RGN@signal_mumuMuMu_withMllWeight@g' combinedMiniFitPlot_temp.C > tempMuMu.C"

dyDatasets=('DYAMC' 'DYPOWHEG')
for q in ${!dyDatasets[*]}
do
	eval "sed 's@MODE@${dyDatasets[$q]}@g' tempEE.C > combinedMiniFitPlot.C"
	eval "root -l -b -q combinedMiniFitPlot.C"
	rm combinedMiniFitPlot.C

	eval "sed 's@MODE@${dyDatasets[$q]}@g' tempMuMu.C > combinedMiniFitPlot.C"
	eval "root -l -b -q combinedMiniFitPlot.C"
	rm combinedMiniFitPlot.C
done
rm temp*.C 


#modify combinedMiniPlotterOneProcessMuMu.C and runCombinedMiniPlotterOneProcessMuMu.C to make plots only with
#ttbar in EE, EMu, and MuMu signal regions, DYAMC in EE and MuMu signal regions
#and real emu data in EMu signal region
#these plots will be used to check the quality of fits done by analysis.cpp
eval "sed 's@MuMu@EMu@g' combinedMiniPlotterOneProcessMuMu.C > tempSignalEMu.C"
eval "sed 's@RGN@flavoursidebandEMu@g' tempSignalEMu.C > tempSignalEMuTwo.C"
#now process tempSignalEMuTwo.C into a form which can be used with real data or TTBar
#real data EMu
eval "sed 's@MODE@data@g' tempSignalEMuTwo.C > tempSignalEMuData.C"
eval "sed 's@OneProcess@RealData@g' tempSignalEMuData.C > combinedMiniPlotterRealDataEMu.C"
eval "sed 's@OneProcessMuMu@RealDataEMu@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterRealDataEMu.C"
#TTBar MC EMu
eval "sed 's@MODE@TT@g' tempSignalEMuTwo.C > tempSignalEMuTT.C"
eval "sed 's@OneProcess@TTBar@g' tempSignalEMuTT.C > combinedMiniPlotterTTBarEMu.C"
eval "sed 's@OneProcessMuMu@TTBarEMu@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterTTBarEMu.C"

#now process combinedMiniPlotterOneProcessMuMu.C into a form which can be used for TTBar and DY MuMu
eval "sed 's@RGN@signal_mumuMuMu_withMllWeight@g' combinedMiniPlotterOneProcessMuMu.C > tempSignalMuMu.C"
#DYAMC MuMu
eval "sed 's@MODE@DYAMC@g' tempSignalMuMu.C > tempSignalDYAMCMuMu.C"
eval "sed 's@OneProcess@DYAMC@g' tempSignalDYAMCMuMu.C > combinedMiniPlotterDYAMCMuMu.C"
eval "sed 's@OneProcessMuMu@DYAMCMuMu@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterDYAMCMuMu.C"

#DYPOWHEG MuMu
eval "sed 's@MODE@DYPOWHEG@g' tempSignalMuMu.C > tempSignalDYPOWHEGMuMu.C"
eval "sed 's@OneProcess@DYPOWHEG@g' tempSignalDYPOWHEGMuMu.C > combinedMiniPlotterDYPOWHEGMuMu.C"
eval "sed 's@OneProcessMuMu@DYPOWHEGMuMu@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterDYPOWHEGMuMu.C"

#TTBar MuMu
eval "sed 's@MODE@TT@g' tempSignalMuMu.C > tempSignalTTMuMu.C"
eval "sed 's@_withMllWeight@@g' tempSignalTTMuMu.C > tempSignalCorrTTMuMu.C"
eval "sed 's@OneProcess@TT@g' tempSignalCorrTTMuMu.C > combinedMiniPlotterTTMuMu.C"
eval "sed 's@OneProcessMuMu@TTMuMu@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterTTMuMu.C"

#now process combinedMiniPlotterOneProcessMuMu.C into a form which can be used for TTBar and DY EE
eval "sed 's@RGN@signal_eeEE_withMllWeight@g' combinedMiniPlotterOneProcessMuMu.C > tempSignalEEPrecursor.C"
eval "sed 's@MuMu@EE@g' tempSignalEEPrecursor.C > tempSignalEE.C"

#DYAMC EE
eval "sed 's@MODE@DYAMC@g' tempSignalEE.C > tempSignalDYAMCEE.C"
eval "sed 's@OneProcess@DYAMC@g' tempSignalDYAMCEE.C > combinedMiniPlotterDYAMCEE.C"
eval "sed 's@OneProcessMuMu@DYAMCEE@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterDYAMCEE.C"

#DYPOWHEG EE
eval "sed 's@MODE@DYPOWHEG@g' tempSignalEE.C > tempSignalDYPOWHEGEE.C"
eval "sed 's@OneProcess@DYPOWHEG@g' tempSignalDYPOWHEGEE.C > combinedMiniPlotterDYPOWHEGEE.C"
eval "sed 's@OneProcessMuMu@DYPOWHEGEE@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterDYPOWHEGEE.C"

#TTBar EE
eval "sed 's@MODE@TT@g' tempSignalEE.C > tempSignalTTEE.C"
eval "sed 's@_withMllWeight@@g' tempSignalTTEE.C > tempSignalCorrTTEE.C"
eval "sed 's@OneProcess@TT@g' tempSignalCorrTTEE.C > combinedMiniPlotterTTEE.C"
eval "sed 's@OneProcessMuMu@TTEE@g' runCombinedMiniPlotterOneProcessMuMu.C > runCombinedMiniPlotterTTEE.C"

eval "root -l -b -q runCombinedMiniPlotterRealDataEMu.C"
eval "root -l -b -q runCombinedMiniPlotterTTBarEMu.C"
eval "root -l -b -q runCombinedMiniPlotterDYAMCMuMu.C"
eval "root -l -b -q runCombinedMiniPlotterDYPOWHEGMuMu.C"
eval "root -l -b -q runCombinedMiniPlotterTTMuMu.C"
eval "root -l -b -q runCombinedMiniPlotterDYAMCEE.C"
eval "root -l -b -q runCombinedMiniPlotterDYPOWHEGEE.C"
eval "root -l -b -q runCombinedMiniPlotterTTEE.C"

rm tempSignal*.C
rm combinedMiniPlotterRealDataEMu.C combinedMiniPlotterTT*.C combinedMiniPlotterDY*.C 
rm runCombinedMiniPlotterRealDataEMu.C runCombinedMiniPlotterTT*.C runCombinedMiniPlotterDY*.C 

eval "cd ../."
echo "now in cmsWR/."


