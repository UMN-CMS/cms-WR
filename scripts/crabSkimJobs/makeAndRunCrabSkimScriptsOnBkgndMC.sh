#!/bin/bash

#all arrays must have the same size 
datasets=(DYJetsMGrtThn50AmcNloReMiniAODToEETwoHEEP25nsNov11 DYJetsMGrtThn50AmcNloReMiniAODToEMuHEEPandIsHighPtID25nsNov11 DYJetsMGrtThn50MadgraphReMiniAODToEETwoHEEP25nsNov11 DYJetsMGrtThn50MadgraphReMiniAODToEMuHEEPandIsHighPtID25nsNov11)
channel=(EEJJ EMuJJ EEJJ EMuJJ)
inputData=('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM')

#datasets=(DYJetsM100To200ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM100To200ReMiniAODToEETwoHEEP25nsOct15)
#channel=(EMuJJ EEJJ)
#inputData=('/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM')

#datasets=(DYJetsM200To400ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM200To400ReMiniAODToEETwoHEEP25nsOct15 DYJetsM400To500ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM400To500ReMiniAODToEETwoHEEP25nsOct15 DYJetsM500To700ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM500To700ReMiniAODToEETwoHEEP25nsOct15 DYJetsM700To800ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM700To800ReMiniAODToEETwoHEEP25nsOct15 DYJetsM800To1000ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM800To1000ReMiniAODToEETwoHEEP25nsOct15 DYJetsM1500To2000ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM1500To2000ReMiniAODToEETwoHEEP25nsOct15 DYJetsM2000To3000ReMiniAODToEMuHEEPandIsHighPtID25nsOct15 DYJetsM2000To3000ReMiniAODToEETwoHEEP25nsOct15 WJetsReMiniAODToEMuHEEPandIsHighPtID25nsOct15 WJetsReMiniAODToEETwoHEEP25nsOct15 TopWReMiniAODToEMuHEEPandIsHighPtID25nsOct15 TopWReMiniAODToEETwoHEEP25nsOct15 AntiTopWReMiniAODToEMuHEEPandIsHighPtID25nsOct15 AntiTopWReMiniAODToEETwoHEEP25nsOct15 TTPowhegPythiaReMiniAODToEMuHEEPandIsHighPtID25nsOct15 TTPowhegPythiaReMiniAODToEETwoHEEP25nsOct15)
#channel=(EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ EMuJJ EEJJ)
#inputData=('/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM' '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM' '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM' '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM' '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1/MINIAODSIM' '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1/MINIAODSIM')

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


