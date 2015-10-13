#!/bin/bash

#bkgndProcess must match an existing folder name in /eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/
#bkgndProcess=ttBar
bkgndProcess=dyPlusJets

#q btwn 0 and 275 for ttbar
#q btwn 289 and 501 for dyPlusJets
for q in {289..501} 
do
	#replace BBKKGGNNDD and NUM in .csh file and bkgndTemp.py file
	#uncomment two input file paths in bkgndTemp.py per loop iteration
	#replace BBKKGGNNDD and NUM in job submission file
	f=$((2*q+92))
	l=$((2*q+1+92))
	
	eval "sed 's/NUM/$q/g' bkgndTemp.py > recoTempOne.py"
	eval "sed 's/BBKKGGNNDD/$bkgndProcess/g' recoTempOne.py > recoTempTwo.py"
	eval "sed '${f},${l}s/#//' recoTempTwo.py > bkgndRecoAnalysis_${bkgndProcess}_$q.py"

	eval "sed 's/NUM/$q/g' runRecoAnalysis.csh > runTempOne.csh"
	eval "sed 's/BBKKGGNNDD/$bkgndProcess/g' runTempOne.csh > runRecoAnalysis_${bkgndProcess}_$q.csh"
	
	eval "sed 's/NUM/$q/g' recoAnalysis > recoAnalysisTemp"
	eval "sed 's/BBKKGGNNDD/$bkgndProcess/g' recoAnalysisTemp > recoAnalysis_${bkgndProcess}_$q"
	
	rm recoTempOne.py recoTempTwo.py runTempOne.csh recoAnalysisTemp

done


