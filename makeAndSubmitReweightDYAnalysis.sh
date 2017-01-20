#!/bin/bash

##AMC
mode='DYAMC'

##MadIncl LowHT
#mode='DYMAD'

##MadHT all bins
#mode='DYMADHT'

for i in {0..100}
do
	#replace WGTNUM by an integer, MODE by a string, and CHNL by a string
	eval "sed 's@WGTNUM@$i@g' runAnalysisReweightDYMC.sh > tempOne.sh"
	eval "sed 's@MODE@$mode@g' tempOne.sh > tempTwo.sh"
	eval "sed 's@CHNL@EE@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
	eval "sed 's@CHNL@MuMu@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
	eval "chmod u+x runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
	eval "chmod u+x runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
	rm tempOne.sh tempTwo.sh

	#submit the job(s)
	#echo "bsub -R 'pool>2000' -q 8nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
	#echo "bsub -R 'pool>2000' -q 8nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
	eval "bsub -R 'pool>2000' -q 8nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
	eval "bsub -R 'pool>2000' -q 8nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
	
done


