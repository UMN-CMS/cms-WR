#!/bin/bash

#provide a list of input and output directories and file names, and this script will rename the files, hadd indvidiual groups
#of files into single files, and copy the hadd'd files to a different directory

#dir, wrMass, and nuMass must have the same number of elements
commonDir='/eos/uscms/store/user/skalafut'

#directory where the hadd'd files should be moved to
copyDir='/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region'


#the elements in dir should not begin or end with a forward slash
#dir=('WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_800_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_093442/0000')
#wrMass=(800)
#nuMass=(400)

dir=('WRToNuEToEEJJ_MW-1400_MNu-700_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_1400_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095455/0000')
wrMass=(1400)
nuMass=(700)

#dir=('WRToNuEToEEJJ_MW-1200_MNu-600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_1200_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095441/0000'  'WRToNuEToEEJJ_MW-1400_MNu-700_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_1400_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095455/0000'  'WRToNuEToEEJJ_MW-1600_MNu-800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_1600_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095510/0000'  'WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_2000_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095546/0000'  'WRToNuEToEEJJ_MW-2400_MNu-1200_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_2400_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095613/0000'  'WRToNuEToEEJJ_MW-2600_MNu-1300_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_2600_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095626/0000'  'WRToNuEToEEJJ_MW-2800_MNu-1400_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_2800_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095644/0000'  'WRToNuEToEEJJ_MW-3000_MNu-1500_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_3000_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095659/0000'  'WRToNuEToEEJJ_MW-3200_MNu-1600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_3200_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095712/0000'  'WRToNuEToEEJJ_MW-3600_MNu-1800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_3600_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095726/0000'  'WRToNuEToEEJJ_MW-3800_MNu-1900_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_3800_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095740/0000'  'WRToNuEToEEJJ_MW-4000_MNu-2000_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_4000_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095758/0000'  'WRToNuEToEEJJ_MW-4200_MNu-2100_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_4200_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095817/0000'  'WRToNuEToEEJJ_MW-4400_MNu-2200_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_4400_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095832/0000'  'WRToNuEToEEJJ_MW-4600_MNu-2300_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_4600_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095846/0000'  'WRToNuEToEEJJ_MW-4800_MNu-2400_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_4800_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095904/0000'  'WRToNuEToEEJJ_MW-5000_MNu-2500_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_5000_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095924/0000'  'WRToNuEToEEJJ_MW-5200_MNu-2600_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_5200_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095938/0000'  'WRToNuEToEEJJ_MW-5600_MNu-2800_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_5600_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_095952/0000'  'WRToNuEToEEJJ_MW-5800_MNu-2900_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_5800_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_100010/0000'  'WRToNuEToEEJJ_MW-6000_MNu-3000_TuneCUETP8M1_13TeV-pythia8/WRToEEJJ_HLTAndRecoOffline_MWR_6000_MNu_halfMWR_25ns_Nov24_13TeV_analyzed_eejj/151124_100028/0000')
#wrMass=(1200 1400 1600 2000 2400 2600 2800 3000 3200 3600 3800 4000 4200 4400 4600 4800 5000 5200 5600 5800 6000)
#nuMass=(600 700 800 1000 1200 1300 1400 1500 1600 1800 1900 2000 2100 2200 2300 2400 2500 2600 2800 2900 3000)


for r in ${!dir[*]}
do
	for q in {1..12}
	do
		#rename the TTree files to display the WR and Nu masses
		mv $commonDir/${dir[$r]}/analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_$q.root $commonDir/${dir[$r]}/analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_$q.root
		#cp $commonDir/${dir[$r]}/analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_$q.root $commonDir/${dir[$r]}/analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_$q.root
	
	done

	#echo "finished renaming one group of files"

	#execute hadd
	eval "hadd -k -O $commonDir/${dir[$r]}/all_analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}.root $commonDir/${dir[$r]}/analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}_*.root"
	#echo "finished hadding one group of files"
	
	#move the hadd'd file to copyDir
	mv $commonDir/${dir[$r]}/all_analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion_WR_M-${wrMass[$r]}_Nu_M-${nuMass[$r]}.root $copyDir/.
	#echo "finished copying one hadd file to copyDir"

done

