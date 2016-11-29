#!/bin/bash
#run this script from cmsWR dir, not scripts dir
mass_n=$(seq 1 25)

eval "make"
eval "mkdir -p plots"

for mass in $mass_n
do
	eval "./bin/analysis -c EE -m signal --makeSelectorPlots true --signalN $mass >& wrSignalEEWithSelHists_$mass.txt"
	eval "./bin/analysis -c MuMu -m signal --makeSelectorPlots true --signalN $mass >& wrSignalMuMuWithSelHists_$mass.txt"
done

##EE dy
#bsub -q $queue -eo $logdir/dy_ee_${tag}.err -oo $logdir/dy_ee_${tag}.out "$PWD/runJob.sh EE DYAMC $finaldir $PWD \"$arguments\""

##MUMU dy
#bsub -q $queue -eo $logdir/dy_mumu_${tag}.err -oo $logdir/dy_mumu_${tag}.out "$PWD/runJob.sh MuMu DYAMC $finaldir $PWD \"$arguments\""

##EE ttbar mc
#bsub -q $queue -eo $logdir/ttbar_ee_${tag}.err -oo $logdir/ttbar_ee_${tag}.out "$PWD/runJob.sh EE TT $finaldir $PWD \"$arguments\""

##MUMU ttbar mc
#bsub -q $queue -eo $logdir/ttbar_mumu_${tag}.err -oo $logdir/ttbar_mumu_${tag}.out "$PWD/runJob.sh MuMu TT $finaldir $PWD \"$arguments\""
