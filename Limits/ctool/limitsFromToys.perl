#!/usr/bin/perl

$grid=$ARGV[0];
$loc=$grid;
$loc=~s|/[^/]+||;
if ($loc eq $grid) {
    $loc=".";
}
$ifile=$grid;
$ifile=~s|.*/||;
$ifile=~s|.root||;
$stub=$ifile;
$ifile="${loc}/input/${ifile}";
$stub=~/([0-9]+)_([0-9]+)/;
$mw=$1;
$mn=$2;
$summary="${loc}/${stub}-summary.log";
open(SUM,">$summary");

print SUM "MW=${mw}\n";
print SUM "MN=${mn}\n";

# now, we extract the limits
limitFromGrid(-1,"obs");
limitFromGrid(0.5,"exp");
limitFromGrid(0.84,"exp_p1s");
limitFromGrid(0.16,"exp_m1s");
limitFromGrid(0.975,"exp_p2s");
limitFromGrid(0.025,"exp_m2s");

close(SUM);

sub limitFromGrid() {
    my($pt,$label)=@_;
    
    if ($pt>0) {
	$cmd="combine -M HybridNew --grid=${grid} --expectedFromGrid ${pt} -m ${mw} ${ifile}";
    } else {
	$cmd="combine -M HybridNew --grid=${grid} -m ${mw} ${ifile}";
    }

    print "===== $cmd \n";
    open(CMD,"$cmd |");
    while (<CMD>) {
	print LOG;
	if (/Limit: r < ([0-9.]+)/) {
	    $rv=$1;
	}
    }
    close(CMD);
    print SUM "${label} = $rv\n";
    
}

