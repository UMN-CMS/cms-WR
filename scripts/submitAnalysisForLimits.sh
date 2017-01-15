arguments="--toys 3200 --ignoreDyScaleFactors false"
logdir=$PWD/log
mass_n=$(seq 1 25)

queue=2nd
tag=3200toysAllSyst
finaldir=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/$tag/

#for mass in $mass_n
#do
#	bsub -q $queue -eo $logdir/${mass}_signal_ee_${tag}.err -oo $logdir/${mass}_signal_ee_${tag}.out "$PWD/runJob.sh EE signal $finaldir $PWD \"$arguments --signalN $mass\""
#	bsub -q $queue -eo $logdir/${mass}_signal_mumu_${tag}.err -oo $logdir/${mass}_signal_mumu_${tag}.out "$PWD/runJob.sh MuMu signal $finaldir $PWD \"$arguments --signalN $mass\""
#done

#EMu data
#bsub -q $queue -eo $logdir/data_emu_ee_${tag}.err   -oo $logdir/data_emu_ee_${tag}.out   "$PWD/runJob.sh EMu data $finaldir $PWD \"$arguments  --cut_channel EE\""
#bsub -q $queue -eo $logdir/data_emu_mumu_${tag}.err -oo $logdir/data_emu_mumu_${tag}.out "$PWD/runJob.sh EMu data $finaldir $PWD \"$arguments  --cut_channel MuMu\""

#DYAMC EE dy
#bsub -q $queue -eo $logdir/dy_ee_${tag}.err -oo $logdir/dy_ee_${tag}.out "$PWD/runJob.sh EE DYAMC $finaldir $PWD \"$arguments\""

#DYAMC MUMU dy
#bsub -q $queue -eo $logdir/dy_mumu_${tag}.err -oo $logdir/dy_mumu_${tag}.out "$PWD/runJob.sh MuMu DYAMC $finaldir $PWD \"$arguments\""

#DYMADHT EE dy
bsub -q $queue -eo $logdir/dy_ee_${tag}.err -oo $logdir/dy_ee_${tag}.out "$PWD/runJob.sh EE DYMADHT $finaldir $PWD \"$arguments\""

#DYMADHT MUMU dy
bsub -q $queue -eo $logdir/dy_mumu_${tag}.err -oo $logdir/dy_mumu_${tag}.out "$PWD/runJob.sh MuMu DYMADHT $finaldir $PWD \"$arguments\""


##EE ttbar mc
#bsub -q $queue -eo $logdir/ttbar_ee_${tag}.err -oo $logdir/ttbar_ee_${tag}.out "$PWD/runJob.sh EE TT $finaldir $PWD \"$arguments\""

##MUMU ttbar mc
#bsub -q $queue -eo $logdir/ttbar_mumu_${tag}.err -oo $logdir/ttbar_mumu_${tag}.out "$PWD/runJob.sh MuMu TT $finaldir $PWD \"$arguments\""
