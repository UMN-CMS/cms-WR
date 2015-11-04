#!/bin/bash

#provide a list of input and output directories and file names, and this script will rename the files, hadd indvidiual groups
#of files into single files, and copy the hadd'd files to a different directory

#dir, wrMass, and nuMass must have the same number of elements
commonDir='/eos/uscms/store/user/skalafut'

#the elements in dir should not begin or end with a forward slash
#dir=('WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_800_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_093933/0000')
#wrMass=(800)
#nuMass=(400)

dir=('WRToNuEToEEJJ_MW-1000_MNu-500_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_1000_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094002/0000'  'WRToNuEToEEJJ_MW-1200_MNu-600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_1200_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094049/0000'  'WRToNuEToEEJJ_MW-1400_MNu-700_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_1400_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094102/0000'  'WRToNuEToEEJJ_MW-1600_MNu-800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_1600_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094115/0000'  'WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_2000_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094143/0000'  'WRToNuEToEEJJ_MW-2200_MNu-1100_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_2200_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094203/0000'  'WRToNuEToEEJJ_MW-2400_MNu-1200_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_2400_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094216/0000'  'WRToNuEToEEJJ_MW-2600_MNu-1300_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_2600_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094230/0000'  'WRToNuEToEEJJ_MW-2800_MNu-1400_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_2800_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094244/0000'  'WRToNuEToEEJJ_MW-3000_MNu-1500_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_3000_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094304/0000'  'WRToNuEToEEJJ_MW-3200_MNu-1600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_3200_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094318/0000'  'WRToNuEToEEJJ_MW-3600_MNu-1800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_3600_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094333/0000'  'WRToNuEToEEJJ_MW-3800_MNu-1900_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_3800_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094346/0000'  'WRToNuEToEEJJ_MW-4400_MNu-2200_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_4400_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094431/0000'  'WRToNuEToEEJJ_MW-5000_MNu-2500_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_5000_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094517/0000'  'WRToNuEToEEJJ_MW-5200_MNu-2600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_5200_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094532/0000'  'WRToNuEToEEJJ_MW-5600_MNu-2800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_5600_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094547/0000'  'WRToNuEToEEJJ_MW-5800_MNu-2900_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_5800_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094600/0000'  'WRToNuEToEEJJ_MW-6000_MNu-3000_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_CheckGenMWR_MWR_6000_MNu_halfMWR_25ns_Nov03_13TeV_analyzed_eejj/151103_094627/0000')
wrMass=(1000 1200 1400 1600 2000 2200 2400 2600 2800 3000 3200 3600 3800 4400 5000 5200 5600 5800 6000)
nuMass=(500 600 700 800 1000 1100 1200 1300 1400 1500 1600 1800 1900 2200 2500 2600 2800 2900 3000)

copyDir='/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass'

for r in ${!dir[*]}
do
	for q in {1..12}
	do
		mv $commonDir/${dir[$r]}/genWrNuAndDecayKinematicsNoMatchingInfo_$q.root $commonDir/${dir[$r]}/genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_$q.root
		#cp $commonDir/${dir[$r]}/genWrNuAndDecayKinematicsNoMatchingInfo_$q.root $commonDir/${dir[$r]}/genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_$q.root
	
	done

	#echo "finished renaming one group of files"

	#execute the hadd call
	eval "hadd -k -O $commonDir/${dir[$r]}/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}.root $commonDir/${dir[$r]}/genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_*.root"
	#echo "finished hadding one group of files"
	
	#move the hadd'd file to a different directory
	cp $commonDir/${dir[$r]}/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}.root $copyDir/.
	#echo "finished copying one hadd file to a different directory"

done

