#!/usr/bin/perl

use Getopt::Long;

$special=0;
$xsec=1.0;

GetOptions(	   "special=i" => \$special);

foreach $file (@ARGV) {

    $file=~/_([0-9]+)_([0-9]+).log/;
    $mw=$1;
    $mn=$2;
    
    open(OF,$file);
    while (<OF>) {

	if (/xsec=([0-9.]+)/) {
	    $xsec=$1;
	}

	$mode=1 if (/Observed/);
	$mode=2 if (/Expected/);
	
	if (/Limit: r < ([0-9.]+) /) {
	    $obs=$1 if ($mode==1);
	}
	if (/median expected limit: r < ([0-9.]+) /) {
	    $exp=$1 if ($mode==2);
	}
	if (/68%/ && /([0-9.]+) < r < ([0-9.]+)/) {
	    $exp_m1s=$1;
	    $exp_p1s=$2;
	}
	if (/95% exp/ && /([0-9.]+) < r < ([0-9.]+)/) {
	    $exp_m2s=$1;
	    $exp_p2s=$2;
	}
    }
    close(OF);
    
    printf("%5d %5d %9.7f %9.7f %9.7f %9.7f %9.7f %9.7f %3d\n",
	   $mw, $mn,
	   $obs*$xsec,
	   $exp*$xsec,$exp_m2s*$xsec,$exp_m1s*$xsec,$exp_p1s*$xsec,$exp_p2s*$xsec,
	   $special);
}
