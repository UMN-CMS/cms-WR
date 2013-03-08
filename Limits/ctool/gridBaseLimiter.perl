#!/usr/bin/perl

$dfile=$ARGV[0];
$ofile=$ARGV[1];
$log=$ARGV[2];
$stub=$dfile;
$stub=~s/.*\///;
$mw=$stub;
$mw=~s/[a-z]+_([0-9]+)_.*/$1/;
$parallelism=12;
$points=24;
$wd="/export/scratch/users/jmmans/${stub}";

print "$wd\n";
system("rm -rf $wd");
system("mkdir -p $wd");
chdir($wd);

# Seed with the ProfileLikelihood
$cmd="combine -v0 -n wr -m $mw -M ProfileLikelihood -t 400 $dfile|";

open(LOG,">${log}");
open(SEED,$cmd);
while (<SEED>) {
    print LOG;
    if (/median expected limit: r < ([0-9.]+)/) {
	$seed_median=$1;
    }

    if (/68% expected band : ([0-9.]+) < r < ([0-9.]+)/) {
	$seed_68l=$1;
	$seed_68h=$2;
    }

    if (/95% expected band : ([0-9.]+) < r < ([0-9.]+)/) {
	$seed_95l=$1;
	$seed_95h=$2;
    }

}
close(SEED);
print LOG "$seed_median $seed_68l $seed_68h\n";

$lp=log($seed_95l*0.75);
$hp=log($seed_95h*1.5);

$running=0;
$files="";
for ($i=0; $i<$points; $i++) {
    $rv=$lp+$i*($hp-$lp)/$points;
    $rv=exp($rv);
    print LOG "$i $rv\n";
    $seed=$i+248100;

    $cmd="combine -v0 -n wr -m $mw -M HybridNew --clsAcc 0 -t 200 --saveToys --saveHybridResult --singlePoint $rv -s $seed --frequentist --testStat LHC $dfile > ${i}.log 2>&1 ";
    $files.=sprintf("higgsCombinewr.HybridNew.mH%d.%d.root ",$mw,$seed);
    if (fork()==0) {
	system($cmd);
	exit(0);
    } 
    $running++;
    while ($running>=$parallelism) {
	
#	print "Waiting\n";
	$av=wait();
	$running-- if ($av!=-1);
	break if ($av==-1);
    }
}

while ($running>0) {
	$av=wait();
	$running-- if ($av!=-1);
	break if ($av==-1);
}

close(LOG);

unlink($ofile);

system("hadd $ofile $files >> $log 2>&1 ");

