#!/bin/bash
#run this script starting in cmsWR/scripts/.
#this script can only be used at LXPLUS, and after a grid proxy has been moved into the lxplus ~/. user area
#this script launches jobs on the lxplus batch system to run analysis.cpp over real data and DY MC dytagandprobe minitrees
#main user modified parameter in this script identifies whether or not to ignore the DY MLL scaling factors
#which are calculated in events with at least two leptons passing ID with pt>35, and at least two jets passing tight ID with pt>40
ignoreDYMLLSF="true"

#these user modified parameters identify the path name to the local directory, and the path to the grid proxy
userProxyPath="/afs/cern.ch/user/s/skalafut/x509up_u38430"
localDirPath="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment/CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR"

eval "sed 's@UPDATE@$ignoreDYMLLSF@' <runTagAndProbeTemp.sh >tempA.sh"
eval "sed 's@PROXYPATH@$userProxyPath@' <tempA.sh >tempB.sh"
eval "sed 's@LOCALPATH@$localDirPath@' <tempB.sh >tempC.sh"

datasets=('data' 'DYPOWHEG' 'DYAMC' 'DYMADHT')
chnls=('EE' 'MuMu')
for r in ${!datasets[*]}
do
	for m in ${!chnls[*]}
	do
		#replace MODE and CHNL in tempC.sh
		eval "sed 's@MODE@${datasets[$r]}@' <tempC.sh >tempD.sh"
		eval "sed 's@CHNL@${chnls[$m]}@' <tempD.sh >runTagAndProbe_${chnls[$m]}_${datasets[$r]}Job.sh"

		#now submit the job to lxplus batch system
		#move up one directory for submission, then move back into scripts directory
		eval "cd ../."
		#echo "bsub -R 'rusage[mem=8000]' -q 8nh -J TnP_${chnls[$m]}_${datasets[$r]} < scripts/runTagAndProbe_${chnls[$m]}_${datasets[$r]}Job.sh"
		eval "bsub -R 'rusage[mem=8000]' -q 8nh -J TnP_${chnls[$m]}_${datasets[$r]} < scripts/runTagAndProbe_${chnls[$m]}_${datasets[$r]}Job.sh"
		eval "cd scripts"

		rm tempD.sh
	done

done

rm tempA.sh tempB.sh tempC.sh
