#!/bin/bash

#provide a list of input and output directories and file names, and this script will rename the files, hadd indvidiual groups
#of files into single files, and copy the hadd'd files to a different directory

#dir, wrMass, and nuMass must have the same number of elements
commonDir='/eos/uscms/store/user/skalafut'

#directory where the hadd'd files should be moved to
copyDir='/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass'

#the elements in dir should not begin or end with a forward slash
#dir and stageOneHaddOutputName must have the same number of elements
#dir=('TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0000'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0001'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0002'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0003'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0004'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0005'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0006'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0007'  'TT_TuneCUETP8M1_13TeV-powheg-pythia8/TTBarPowhegPythiaToEEJJ_GenOffline_25ns_Nov19_13TeV_analyzed_eejj/151120_085735/0008')
#stageOneHaddOutputName=('all_genTTBarDecayKinematics_0.root'  'all_genTTBarDecayKinematics_1.root'  'all_genTTBarDecayKinematics_2.root'  'all_genTTBarDecayKinematics_3.root'  'all_genTTBarDecayKinematics_4.root'  'all_genTTBarDecayKinematics_5.root'  'all_genTTBarDecayKinematics_6.root'  'all_genTTBarDecayKinematics_7.root'  'all_genTTBarDecayKinematics_8.root')
dir=('DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_200To400_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145108/0000'  'DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_400To500_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145121/0000'  'DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_500To700_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145134/0000'  'DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_700To800_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145148/0000'  'DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_800To1000_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145201/0000'  'DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_1500To2000_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145219/0000'  'DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsAMCNLO_M_2000To3000_EEJJ_GenOffline_25ns_Nov20Redo_13TeV_analyzed_eejj/151120_145233/0000')
stageOneHaddOutputName=('all_genDYJetsAMCNLO_M200to400_DecayKinematics.root'  'all_genDYJetsAMCNLO_M400to500_DecayKinematics.root'  'all_genDYJetsAMCNLO_M500to700_DecayKinematics.root'  'all_genDYJetsAMCNLO_M700to800_DecayKinematics.root'  'all_genDYJetsAMCNLO_M800to1000_DecayKinematics.root'  'all_genDYJetsAMCNLO_M1500to2000_DecayKinematics.root'  'all_genDYJetsAMCNLO_M2000to3000_DecayKinematics.root')


#noNumsStageOneHaddOutputName and stageTwoHaddOutputName must have the same number of elements
#noNumsStageOneHaddOutputName=('all_genTTBarDecayKinematics')
#stageTwoHaddOutputName=('complete_genTTBarDecayKinematics.root')



for r in ${!dir[*]}
do
	#the bkgnd files do not need to be renamed
	#for q in {1..12}
	#do
	#	#rename the TTree files to display the WR and Nu masses
	#	mv $commonDir/${dir[$r]}/analyzed_tree_allGenAndRecoOfflineCuts_eejjSignalRegion_$q.root $commonDir/${dir[$r]}/analyzed_tree_allGenAndRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_$q.root
	#
	#done

	#execute hadd
	eval "hadd -k -O $commonDir/${dir[$r]}/${stageOneHaddOutputName[$r]} $commonDir/${dir[$r]}/*.root"
	#echo "hadd -k -O $commonDir/${dir[$r]}/${stageOneHaddOutputName[$r]} $commonDir/${dir[$r]}/*.root"
	
	#move the hadd'd file to copyDir
	mv $commonDir/${dir[$r]}/${stageOneHaddOutputName[$r]} $copyDir/.
	#echo "$commonDir/${dir[$r]}/${stageOneHaddOutputName[$r]} $copyDir/."


done

#hadd the hadd'd files together in copyDir, and after hadd remove the files which were passed as input to hadd
#for s in ${!stageTwoHaddOutputName[*]}
#do
#	eval "hadd -k -O $copyDir/${stageTwoHaddOutputName[$s]} $copyDir/${noNumsStageOneHaddOutputName[$s]}_*.root"
#	eval "rm $copyDir/${noNumsStageOneHaddOutputName[$s]}_*.root"
#
#done


