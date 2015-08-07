#!/bin/bash

#all arrays should have the same size 
#datasets=(ZZToEEJJ50nsSkim ZZToEMuJJ50nsSkim WZToEEJJ50nsSkim WZToEMuJJ50nsSkim WJetsToEEJJ50nsSkim WJetsToEMuJJ50nsSkim)
#channel=(EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ)
#inputData=('/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM')

datasets=(ZZToEMuJJ50nsSkim WZToEEJJ50nsSkim WZToEMuJJ50nsSkim WJetsToEEJJ50nsSkim WJetsToEMuJJ50nsSkim)
channel=(EMuJJ EEJJ EMuJJ EEJJ EMuJJ)
inputData=('/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM')





for q in ${!datasets[*]}
do
	
		#replace FNLST with an element from datasets, CHNL with an element from channel, and INPTDT with an element from inputData
		#in skim_bkgndMC_chnl_temp.py
		eval "sed 's/FNLST/${datasets[$q]}/g' skim_bkgndMC_chnl_temp.py > skim_bkgndMC_chnl_one.py"
		eval "sed 's@INPTDT@${inputData[$q]}@g' skim_bkgndMC_chnl_one.py > skim_bkgndMC_chnl_two.py"
		eval "sed 's/CHNL/${channel[$q]}/g' skim_bkgndMC_chnl_two.py > skim_bkgndMC_signalAndLowMassRegions_${datasets[$q]}.py"
		rm skim_bkgndMC_chnl_one.py skim_bkgndMC_chnl_two.py
	
		#submit jobs to crab using the newly created crab skim .py file
		#echo "crab submit -c skim_bkgndMC_signalAndLowMassRegions_${datasets[$q]}.py"
		eval "crab submit -c skim_bkgndMC_signalAndLowMassRegions_${datasets[$q]}.py"


done


