#!/bin/bash
##ONLY run this script in cmsWR/ directory, not a subdirectory
workingDir=$PWD

##AMC
#mode='DYAMC'

##MadIncl LowHT and MadHT all bins
#mode='DYMadInclAndHT'

##for DYMAD and DYAMC
#for i in {0..100}
#do
#	#replace WGTNUM by an integer, MODE by a string, and CHNL by a string
#	eval "sed 's@WGTNUM@$i@g' runAnalysisReweightDYMC.sh > tempOne.sh"
#	eval "sed 's@WORKDIR@$workingDir@g' tempOne.sh > tempOnePfive.sh"
#	eval "sed 's@MODE@$mode@g' tempOnePfive.sh > tempTwo.sh"
#	eval "sed 's@CHNL@EE@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
#	eval "sed 's@CHNL@MuMu@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
#	eval "chmod u+x runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
#	eval "chmod u+x runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
#	rm tempOne.sh tempOnePfive.sh tempTwo.sh
#
#	#submit the job(s)
#	#echo "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
#	#echo "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
#	eval "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
#	eval "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
#	
#done


##WR signal
mode='signal'
mass_n=$(seq 1 25)

#for WR signal
for j in $mass_n
do

	for i in {0..100}
	do
		#replace WGTNUM by an integer, MODE by a string, and CHNL by a string
		eval "sed 's@WGTNUM@$i@g' runAnalysisReweightDYMC.sh > tempOne.sh"
		eval "sed 's@WORKDIR@$workingDir@g' tempOne.sh > tempOnePfive.sh"
		eval "sed 's@MODE@$mode@g' tempOnePfive.sh > tempTwo.sh"
		eval "sed 's@CHNL@EE@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
		eval "sed 's@CHNL@MuMu@g' tempTwo.sh > runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
		eval "chmod u+x runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
		eval "chmod u+x runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
		rm tempOne.sh tempOnePfive.sh tempTwo.sh

		#submit the job(s)
		#echo "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
		#echo "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"
		eval "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_EE_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_EE_pdfWeight_${i}.sh"
		eval "bsub -R 'pool>1500' -q 1nh -J reweightDYMC_${mode}_MuMu_pdfWeight_${i} < runAnalysisReweightDYMC_${mode}_MuMu_pdfWeight_${i}.sh"

	done

done

