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

		#move to test/ dir and make a skim cfg python file
		#eval "cd ../../test/"
		#eval "sed 's@MASSWR@${WrMasses[$h]}@g' temp_skims_cfg.py > skimsOne.py"
		#eval "sed 's@MASSNU@${NuMasses[$h]}@g' skimsOne.py > skimsTwo.py"
		#eval "sed 's@CHNL@EEJJ@g' skimsTwo.py > skimsThree.py"
		#eval "sed 's@NNN@$jobId@g' skimsThree.py > skimsFour.py"
		#eval "sed 's@TAGNAME@WRtoEEJJ_${WrMasses[$h]}_${NuMasses[$h]}_PrivReco@g' skimsFour.py > skims_cfg_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.py"
		#rm skimsOne.py skimsTwo.py skimsThree.py skimsFour.py
	
		##move back to starting dir
		#eval "cd -"

		#submit job to queue
		#echo "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAod_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.sh"
		eval "bsub -R 'rusage[mem=1500] && pool>2000' -q 1nh -J genSimToHltToRecoToMiniAod_${jobId}_${WrMasses[$h]}_${NuMasses[$h]} < runBatchJob_WR_M-${WrMasses[$h]}_ToLNu_M-${NuMasses[$h]}_EEJJ_${jobId}.sh"
	
		#increment jobId number
		let jobId=jobId+1
	done


done

