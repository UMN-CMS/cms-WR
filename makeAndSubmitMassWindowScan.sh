#!/bin/bash
##ONLY run this script in cmsWR/ directory, not a subdirectory
workingDir=$PWD

mass=(800 1000 1200 1400 1600 1800 2000 2400 2600 2800 3000 3200 3600 3800 4000 4200 4400 4600 4800 5000 5200 5600 5800 6000)
#mass=(2000)


for i in ${!mass[*]}
do
	#replace WORKDIR with an abs path string, WRMASS by a number, and CHNL by a string
	#use ee or mumu for channel
	eval "sed 's@WORKDIR@$workingDir@g' runMassWindowScan.sh > tempOnePfive.sh"
	eval "sed 's@WRMASS@${mass[$i]}@g' tempOnePfive.sh > tempTwo.sh"
	eval "sed 's@CHNL@ee@g' tempTwo.sh > runMassWindowScan_MWR_${mass[$i]}_ee.sh"
	eval "sed 's@CHNL@mumu@g' tempTwo.sh > runMassWindowScan_MWR_${mass[$i]}_mumu.sh"
	eval "chmod u+x runMassWindowScan_MWR_${mass[$i]}_ee.sh"
	eval "chmod u+x runMassWindowScan_MWR_${mass[$i]}_mumu.sh"
	rm tempOnePfive.sh tempTwo.sh

	#replace WRMASS with a number and CHNL with ee or mumu
	eval "cd scripts"
	eval "sed 's@WRMASS@${mass[$i]}@g' limit_window_scan.py > tempOne.sh"
	eval "sed 's@CHNL@ee@g' tempOne.sh > limit_window_scan_MWR_${mass[$i]}_ee.py"
	eval "sed 's@CHNL@mumu@g' tempOne.sh > limit_window_scan_MWR_${mass[$i]}_mumu.py"
	rm tempOne.sh
	eval "cd .."
	
	#submit the job(s)
	#echo "bsub -R 'pool>1500' -q 1nh -J scanMassWindows_ee_MWR_${mass[$i]} < runMassWindowScan_MWR_${mass[$i]}_ee.sh"
	#echo "bsub -R 'pool>1500' -q 1nh -J scanMassWindows_mumu_MWR_${mass[$i]} < runMassWindowScan_MWR_${mass[$i]}_mumu.sh"
	eval "bsub -R 'pool>1500' -q 1nh -J scanMassWindows_ee_MWR_${mass[$i]} < runMassWindowScan_MWR_${mass[$i]}_ee.sh"
	eval "bsub -R 'pool>1500' -q 1nh -J scanMassWindows_mumu_MWR_${mass[$i]} < runMassWindowScan_MWR_${mass[$i]}_mumu.sh"
	
done


