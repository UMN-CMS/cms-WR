arguments="--toys 4000"
logdir=$PWD/log
mass_n=$(seq 1 25)

queue=1nd
tag=4000toys
finaldir=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/$tag/

for mass in $mass_n
do
	bsub -q $queue -eo $logdir/${mass}_signal_ee_${tag}.err -oo $logdir/${mass}_signal_ee_${tag}.out "$PWD/runJob.sh EE signal $finaldir $PWD \"$arguments --signalN $mass\""
	bsub -q $queue -eo $logdir/${mass}_signal_mumu_${tag}.err -oo $logdir/${mass}_signal_mumu_${tag}.out "$PWD/runJob.sh MuMu signal $finaldir $PWD \"$arguments --signalN $mass\""
done

#EMu data
bsub -q $queue -eo $logdir/data_emu_ee_${tag}.err   -oo $logdir/data_emu_ee_${tag}.out   "$PWD/runJob.sh EMu data $finaldir $PWD \"$arguments  --cut_channel EE\""
bsub -q $queue -eo $logdir/data_emu_mumu_${tag}.err -oo $logdir/data_emu_mumu_${tag}.out "$PWD/runJob.sh EMu data $finaldir $PWD \"$arguments  --cut_channel MuMu\""

#EE dy
bsub -q $queue -eo $logdir/dy_ee_${tag}.err -oo $logdir/dy_ee_${tag}.out "$PWD/runJob.sh EE DYAMC $finaldir $PWD \"$arguments\""

#MUMU dy
bsub -q $queue -eo $logdir/dy_mumu_${tag}.err -oo $logdir/dy_mumu_${tag}.out "$PWD/runJob.sh MuMu DYAMC $finaldir $PWD \"$arguments\""

##EE ttbar mc
#bsub -q $queue -eo $logdir/ttbar_ee_${tag}.err -oo $logdir/ttbar_ee_${tag}.out "$PWD/runJob.sh EE TT $finaldir $PWD \"$arguments\""

##MUMU ttbar mc
#bsub -q $queue -eo $logdir/ttbar_mumu_${tag}.err -oo $logdir/ttbar_mumu_${tag}.out "$PWD/runJob.sh MuMu TT $finaldir $PWD \"$arguments\""
