arguments="--makeSelectorPlots"
logdir=$PWD/log
mass_n="1000_500 1200_600 1400_700 1600_800 1800_1400 2000_1000 2200_1100 2400_1200 2600_1300 2800_1400 3000_1500 3200_1600 3600_1800 3800_1900 4000_2000 4200_2100 4400_2200 4600_2300 4800_2400 5000_2500 5200_2600 5600_2800 5800_2900 6000_3000 800_400"

queue=1nd
tag=noTrigger
inputdir=/afs/cern.ch/work/p/phansen/public/wr/noTrigger_minitrees/
finaldir=/afs/cern.ch/work/p/phansen/public/wr/$tag/

for mass in $mass_n
do
	mode=WRtoEEJJ_${mass}
	bsub -q $queue -eo $logdir/${mass}_signal_ee_${tag}.err -oo $logdir/${mass}_signal_ee_${tag}.out "$PWD/runJob.sh EE ${mode} $finaldir $PWD \"$arguments --input ${inputdir}/${mode}_1.root\""

	mode=WRtoMuMuJJ_${mass}
	bsub -q $queue -eo $logdir/${mass}_signal_mumu_${tag}.err -oo $logdir/${mass}_signal_mumu_${tag}.out "$PWD/runJob.sh MuMu ${mode} $finaldir $PWD \"$arguments --input ${inputdir}/${mode}_1.root\""
done


