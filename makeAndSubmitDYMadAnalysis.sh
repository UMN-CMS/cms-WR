#!/bin/bash
#all input files to process
#100to200 inputHTBinnedFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_10_2_kP2.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_11_1_pDt.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_12_1_XBz.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_13_3_911.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_14_2_T4l.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_15_1_pVb.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_16_1_l9I.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_1_2_TRt.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_2_1_KOp.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_3_1_3BN.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_4_2_z5k.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_5_1_3G8.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_6_1_95A.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_7_1_WNZ.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_8_1_Vh3.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_9_1_aQu.root')
#200to400 inputHTBinnedFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_1_1_TCb.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_2_1_2qZ.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_3_1_Dby.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_4_1_UVm.root')
#400to600 inputHTBinnedFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_1_1_NYr.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_2_1_263.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_3_1_vfG.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_4_1_Mtb.root')
#600toInf inputHTBinnedFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_1_1_d7e.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_2_1_2RV.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_3_1_5Uz.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_4_1_oeT.root')


#number used to distinguish different jobs
startingCount=1

##for DYMadHTBinned
#htBin='600toInf'
#for i in ${!inputHTBinnedFiles[*]}
#do
#	#replace NNN by a number and INPTFILE with a file path name from inputHTBinnedFiles
#	eval "cd test"
#	eval "sed 's@NNN@$startingCount@g' temp_genAndRecoDYJets_onlyMadHTBinned_cfg.py > tempOne.py"
#	eval "sed 's@XXX@$htBin@g' tempOne.py > tempTwo.py"
#	eval "sed 's@INPTFILE@${inputHTBinnedFiles[$i]}@g' tempTwo.py > genAndRecoDYJets_onlyMadHTBinned${htBin}_cfg_${startingCount}.py"
#	rm tempOne.py tempTwo.py
#	eval "cd .."
#
#	#replace NNN by a number and xxx (but all caps) by a string representing the HT bin, like 100to200
#	eval "sed 's@XXX@$htBin@g' runAnalysisDYMadHTBinnedMC.sh > temp_runAnalysisDYMadHTBinnedMC.sh"
#	eval "sed 's@NNN@$startingCount@g' temp_runAnalysisDYMadHTBinnedMC.sh > runAnalysisDYMadHTBinnedMC_${htBin}_${startingCount}.sh"
#	eval "chmod u+x runAnalysisDYMadHTBinnedMC_${htBin}_${startingCount}.sh"
#	rm temp_runAnalysisDYMadHTBinnedMC.sh
#
#	#submit the job to 1nh queue and request at least 2 GB of disk
#	eval "bsub -R 'pool>2000' -q 1nh -J runAnalysisDYMadHTBinnedMC_${htBin}_Part_${startingCount} < runAnalysisDYMadHTBinnedMC_${htBin}_${startingCount}.sh"
#	
#	let startingCount=startingCount+1
#done

##end for DYMadHTBinned

##for DYMadIncl
inputInclFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_10_2_KFJ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_11_2_PR1.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_12_2_DuM.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_13_2_yHC.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_14_2_bWL.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_15_2_rK5.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_16_2_yJp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_17_2_Gff.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_18_2_KEW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_19_2_3FO.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_1_2_7zA.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_20_2_hc0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_21_2_nur.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_22_1_SSm.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_2_2_89c.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_3_2_yNe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_4_2_aL0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_5_2_RDU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_6_2_wPy.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_7_2_YMG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_8_2_xNp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_9_2_EMH.root')

for i in ${!inputInclFiles[*]}
do
	#replace NNN by a number and INPTFILE with a file path name from inputInclFiles
	eval "cd test"
	eval "sed 's@NNN@$startingCount@g' temp_genAndRecoDYJets_onlyMadIncl_cfg.py > tempOne.py"
	eval "sed 's@INPTFILE@${inputInclFiles[$i]}@g' tempOne.py > genAndRecoDYJets_onlyMadIncl_cfg_${startingCount}.py"
	rm tempOne.py
	eval "cd .."

	#replace NNN by a number
	eval "sed 's@NNN@$startingCount@g' runAnalysisDYMadInclMC.sh > runAnalysisDYMadInclMC_${startingCount}.sh"
	eval "chmod u+x runAnalysisDYMadInclMC_${startingCount}.sh"

	#submit the job to 1nh queue and request at least 2 GB of disk
	eval "bsub -R 'pool>2000' -q 1nh -J runAnalysisDYMadInclMC_Part_${startingCount} < runAnalysisDYMadInclMC_${startingCount}.sh"
	
	let startingCount=startingCount+1
done

##end for DYMadIncl


