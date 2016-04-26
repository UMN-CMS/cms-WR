#!/bin/bash

#cd to cmsWR/. and run the python script test/miniTreeDumpCheck.py to check the number of entries in minitrees
eval "echo commitNumber >& configs/miniTreeEntries.dat"
eval "git log -1 | grep commit | cut -d ' ' -f 2 >> configs/miniTreeEntries.dat"
eval 'python test/miniTreeDumpCheck.py >> configs/miniTreeEntries.dat'

#now make and run analysis.cpp
#make lib and outputFromAnalysis directories if they don't already exist
eval "mkdir -p lib"
eval "mkdir -p outputFromAnalysis"
eval "make"
tAndPdatasetNames=('data' 'DYMAD' 'DYPOWHEG' 'DYAMC')
for r in ${!tAndPdatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkMuMuDYTandP_${tAndPdatasetNames[$r]}.txt &"
	eval "./bin/analysis --c EE --isTagAndProbe true --m ${tAndPdatasetNames[$r]} >& outputFromAnalysis/checkEEDYTandP_${tAndPdatasetNames[$r]}.txt &"
done

#W dataset is WJets
sidebandDatasetNames=('data' 'DYMAD' 'TT' 'W' 'WZ' 'ZZ')
for q in ${sidebandDatasetNames[*]}
do
	eval "./bin/analysis --c MuMu --isLowDiLepton true --m ${sidebandDatasetNames[$q]} >& outputFromAnalysis/checkMuMuLowDileptonMass_${sidebandDatasetNames[$q]}.txt &"
	eval "./bin/analysis --c EE --isLowDiLepton true --m ${sidebandDatasetNames[$q]} >& outputFromAnalysis/checkEELowDileptonMass_${sidebandDatasetNames[$q]}.txt &"
	eval "./bin/analysis --c EMu --m ${sidebandDatasetNames[$q]} >& outputFromAnalysis/checkEMu_${sidebandDatasetNames[$q]}.txt &"
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

#now run the plotting scripts




