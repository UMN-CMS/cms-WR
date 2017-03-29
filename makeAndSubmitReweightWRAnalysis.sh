#!/bin/bash
##ONLY run this script in cmsWR/ directory, not a subdirectory
workingDir=$PWD

##WR signal
mode='signal'
#mass_n=$(seq 1 25)  #1 is MWR 1 TeV, 25 is MWR 800 GeV
mass_n=$(seq 1 16)


for j in $mass_n
do

	#only minitrees with WR mass up to 4.2 TeV have 200 different pdf and fact/renorm weights per evt
	for i in {0..200}
	do
		#replace WGTNUM by an integer, MODE by a string, WORKDIR by abs path, MMAASS by an integer, and CHNL by a string
		eval "sed 's@WGTNUM@$i@g' runAnalysisReweightWRMC.sh > tempOne.sh"
		eval "sed 's@MMAASS@$j@g' tempOne.sh > tempM.sh"
		eval "sed 's@WORKDIR@$workingDir@g' tempM.sh > tempOnePfive.sh"
		eval "sed 's@MODE@$mode@g' tempOnePfive.sh > tempTwo.sh"
		eval "sed 's@CHNL@EE@g' tempTwo.sh > runAnalysisReweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i}.sh"
		eval "sed 's@CHNL@MuMu@g' tempTwo.sh > runAnalysisReweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i}.sh"
		eval "chmod u+x runAnalysisReweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i}.sh"
		eval "chmod u+x runAnalysisReweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i}.sh"
		rm tempOne.sh tempOnePfive.sh tempTwo.sh tempM.sh

		#submit the job(s)
		#echo "bsub -R 'pool>1100' -q 1nh -J reweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i} < runAnalysisReweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i}.sh"
		#echo "bsub -R 'pool>1100' -q 1nh -J reweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i} < runAnalysisReweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i}.sh"
		eval "bsub -R 'pool>1100' -q 1nh -J reweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i} < runAnalysisReweightWRMC_${mode}_signalMass_${j}_EE_pdfWeight_${i}.sh"
		eval "bsub -R 'pool>1100' -q 1nh -J reweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i} < runAnalysisReweightWRMC_${mode}_signalMass_${j}_MuMu_pdfWeight_${i}.sh"

	done

done

