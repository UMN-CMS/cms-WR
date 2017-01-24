#!/bin/bash

#path names to job directories
jobDirs=('crab/skim/crab_skim_QCD_emenr_pt120to170_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt15to20_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt170to300_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt20to30_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt300toInf_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt30to50_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt50to80_SHv12' 'crab/skim/crab_skim_QCD_emenr_pt80to120_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt1000toInf_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt120to170_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt15to20_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt170to300_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt20to30_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt300to470_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt30to50_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt470to600_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt50to80_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt600to800_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt800to1000_SHv12' 'crab/skim/crab_skim_QCD_muenr_pt80to120_SHv12')
subDirs=('crab_QCD_emenr_pt120to170' 'crab_QCD_emenr_pt15to20' 'crab_QCD_emenr_pt170to300' 'crab_QCD_emenr_pt20to30' 'crab_QCD_emenr_pt300toInf' 'crab_QCD_emenr_pt30to50' 'crab_QCD_emenr_pt50to80' 'crab_QCD_emenr_pt80to120' 'crab_QCD_muenr_pt1000toInf' 'crab_QCD_muenr_pt120to170' 'crab_QCD_muenr_pt15to20' 'crab_QCD_muenr_pt170to300' 'crab_QCD_muenr_pt20to30' 'crab_QCD_muenr_pt300to470' 'crab_QCD_muenr_pt30to50' 'crab_QCD_muenr_pt470to600' 'crab_QCD_muenr_pt50to80' 'crab_QCD_muenr_pt600to800' 'crab_QCD_muenr_pt800to1000' 'crab_QCD_muenr_pt80to120')
space='  '

for i in ${!jobDirs[*]}
do
	eval "crab status -d ${jobDirs[$i]}/${subDirs[$i]}"
	#echo "crab status -d ${jobDirs[$i]}/${subDirs[$i]}"
	eval "echo $space"
done

