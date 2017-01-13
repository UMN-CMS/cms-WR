#!/bin/bash
#all input files to process
inputFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_10_2_KFJ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_11_2_PR1.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_12_2_DuM.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_13_2_yHC.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_14_2_bWL.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_15_2_rK5.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_16_2_yJp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_17_2_Gff.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_18_2_KEW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_19_2_3FO.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_1_2_7zA.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_20_2_hc0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_21_2_nur.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_22_1_SSm.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_2_2_89c.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_3_2_yNe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_4_2_aL0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_5_2_RDU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_6_2_wPy.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_7_2_YMG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_8_2_xNp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_9_2_EMH.root')

#number used to distinguish different jobs
startingCount=1

for i in ${!inputFiles[*]}
do
	#replace NNN by a number and INPTFILE with a file path name from inputFiles in temp_runAnalysis_cfg.py
	eval "cd test"
	eval "sed 's@NNN@$startingCount@g' temp_runAnalysis_cfg.py > tempOne.py"
	eval "sed 's@INPTFILE@${inputFiles[$i]}@g' tempOne.py > runAnalysis_cfg_${startingCount}.py"
	rm tempOne.py
	eval "cd .."

	#replace NNN by a number in runAnalysisDYMC.sh
	eval "sed 's@NNN@$startingCount@g' runAnalysisDYMC.sh > runAnalysisDYMC_${startingCount}.sh"
	eval "chmod u+x runAnalysisDYMC_${startingCount}.sh"

	#submit the job to 1nh queue and request at least 2 GB of disk
	eval "bsub -R 'pool>2000' -q 1nh -J runAnalysisDYMCPart_${startingCount} < runAnalysisDYMC_${startingCount}.sh"
	
	let startingCount=startingCount+1
done

