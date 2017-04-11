#!/bin/bash
#execute this script in its local directory

#WrMasses=(2400 2400 4000 4000)
#NuMasses=(2200 500 3800 800)
WrMasses=(2400)
NuMasses=(2200)
jobStartingDir="$PWD"
gridProxyPath="/afs/cern.ch/user/s/skalafut/x509up_u38430"
nJobs=1	#each job produces 250 events
jobId=1	#increment by 1 after submitting each job

for h in ${!WrMasses[*]}
do
	while [ $jobId -le $nJobs ]
	do

		#############EEJJ#################3
		#make python config file for GEN-SIM step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_GEN-SIM_EEJJ.py > WR_tempEEJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempEEJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempEEJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN-SIM_EEJJ_${jobId}.py"
		rm WR_tempEEJJ.py

		#make python config file for HLT step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_HLT_EEJJ.py > WR_tempEEJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempEEJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempEEJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_EEJJ_${jobId}.py"
		rm WR_tempEEJJ.py

		#make python config file for RECO step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_RECO_EEJJ.py > WR_tempEEJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempEEJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempEEJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_RECO_EEJJ_${jobId}.py"
		rm WR_tempEEJJ.py

		#make python config file for MINIAOD step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD_EEJJ.py > WR_tempEEJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempEEJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempEEJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_EEJJ_${jobId}.py"
		rm WR_tempEEJJ.py


		#make bash script which tells batch job what to do
		eval "sed 's@LOCALPATH@$jobStartingDir@' runBatchJobTempEEJJ.sh > jobOneEEJJ.sh"
		eval "sed 's@PROXYPATH@$gridProxyPath@' jobOneEEJJ.sh > jobTwoEEJJ.sh"
		eval "sed -i 's@NUM@$jobId@g' jobTwoEEJJ.sh"
		eval "sed 's@MASSWR@${WrMasses[$h]}@g' jobTwoEEJJ.sh > jobThreeEEJJ.sh"
		eval "sed 's@MASSNU@${NuMasses[$h]}@g' jobThreeEEJJ.sh > runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.sh"
		rm jobOneEEJJ.sh jobTwoEEJJ.sh jobThreeEEJJ.sh

		#submit job to queue
		#echo "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAodToSkimToMinitree_EEJJ_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.sh"
		eval "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAodToSkimToMinitree_EEJJ_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.sh"


		#############MuMuJJ#################3
		#make python config file for GEN-SIM step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_GEN-SIM_MuMuJJ.py > WR_tempMuMuJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempMuMuJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempMuMuJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_GEN-SIM_MuMuJJ_${jobId}.py"
		rm WR_tempMuMuJJ.py

		#make python config file for HLT step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_HLT_MuMuJJ.py > WR_tempMuMuJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempMuMuJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempMuMuJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_HLT_MuMuJJ_${jobId}.py"
		rm WR_tempMuMuJJ.py

		#make python config file for RECO step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_RECO_MuMuJJ.py > WR_tempMuMuJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempMuMuJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempMuMuJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_RECO_MuMuJJ_${jobId}.py"
		rm WR_tempMuMuJJ.py

		#make python config file for MINIAOD step                     
		eval "sed 's/MASSWR/${WrMasses[$h]}/g'  WR_M-UNDEF_ToLNu_M-UNDEF_miniAOD_MuMuJJ.py > WR_tempMuMuJJ.py"   
		eval "sed -i 's@NUM@$jobId@g' WR_tempMuMuJJ.py"
		eval "sed 's/MASSNU/${NuMasses[$h]}/g'  WR_tempMuMuJJ.py > WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_miniAOD_MuMuJJ_${jobId}.py"
		rm WR_tempMuMuJJ.py


		#make bash script which tells batch job what to do
		eval "sed 's@LOCALPATH@$jobStartingDir@' runBatchJobTempMuMuJJ.sh > jobOneMuMuJJ.sh"
		eval "sed 's@PROXYPATH@$gridProxyPath@' jobOneMuMuJJ.sh > jobTwoMuMuJJ.sh"
		eval "sed -i 's@NUM@$jobId@g' jobTwoMuMuJJ.sh"
		eval "sed 's@MASSWR@${WrMasses[$h]}@g' jobTwoMuMuJJ.sh > jobThreeMuMuJJ.sh"
		eval "sed 's@MASSNU@${NuMasses[$h]}@g' jobThreeMuMuJJ.sh > runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_MuMuJJ_${jobId}.sh"
		rm jobOneMuMuJJ.sh jobTwoMuMuJJ.sh jobThreeMuMuJJ.sh

		#submit job to queue
		#echo "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAodToSkimToMinitree_MuMuJJ_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_MuMuJJ_${jobId}.sh"
		eval "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAodToSkimToMinitree_MuMuJJ_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_MuMuJJ_${jobId}.sh"


		#increment jobId number
		let jobId=jobId+1
	done


done

