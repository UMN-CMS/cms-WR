#!/bin/bash

#all arrays must have the same size 
datasets=(TTBarToEMuHEEPandIsHighPtID25nsSept7 DYJetsMadgraphToEMuHEEPandIsHighPtID25nsSept7 ZZToEMuHEEPandIsHighPtID25nsSept7 WZToEMuHEEPandIsHighPtID25nsSept7 WJetsToEMuHEEPandIsHighPtID25nsSept7 TTBarToEETwoHEEP25nsSept7 DYJetsMadgraphToEETwoHEEP25nsSept7 ZZToEETwoHEEP25nsSept7 WZToEETwoHEEP25nsSept7 WJetsToEETwoHEEP25nsSept7)
channel=(EMuJJ EMuJJ EMuJJ EMuJJ EMuJJ EEJJ EEJJ EEJJ EEJJ EEJJ)
inputData=('/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1/MINIAODSIM'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'  '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'  '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'  '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1/MINIAODSIM'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'  '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'  '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'  '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM')


for q in ${!datasets[*]}
do
	
		#replace FNLST with an element from datasets, CHNL with an element from channel, and INPTDT with an element from inputData
		#in skim_bkgndMC_chnl_temp.py
		eval "sed 's/FNLST/${datasets[$q]}/g' skim_bkgndMC_chnl_temp.py > skim_bkgndMC_chnl_one.py"
		eval "sed 's@INPTDT@${inputData[$q]}@g' skim_bkgndMC_chnl_one.py > skim_bkgndMC_chnl_two.py"
		eval "sed 's/CHNL/${channel[$q]}/g' skim_bkgndMC_chnl_two.py > skim_bkgndMC_${datasets[$q]}.py"
		rm skim_bkgndMC_chnl_one.py skim_bkgndMC_chnl_two.py
	
		#submit jobs to crab using the newly created crab skim .py file
		#echo "crab submit -c skim_bkgndMC_${datasets[$q]}.py"
		eval "crab submit -c skim_bkgndMC_${datasets[$q]}.py"


done


