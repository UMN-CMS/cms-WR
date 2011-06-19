#!/usr/bin/perl

use Getopt::Long;

$toys=100;
$method="MarkovChainMC";
$jobBase="default";


GetOptions("toys=i" => \$toys,
	   "method=s" => \$method,
	   "jobname=s" => \$jobBase);



$pwd=`pwd`;
chomp($pwd);
$exe=$pwd."/scanPoint.sh";
$workloc="/local/cms/user/".$ENV{"USER"}."/heavyNuLimits";
if ($jobBase ne "default") {
    $workloc=$workloc."/${jobBase}";
}

$zone2010="/local/cms/user/jmmans/heavyNu_jun14/2010";

$zone2011="/local/cms/user/jmmans/heavyNu_jun14/2011";

$lumi2011=204;

$fileLoc=$workloc."/input";
system("mkdir -p $fileLoc");
system("rm -rf ${workloc}/roostats*");

# We use the 2011 as the master list of points
open(PLIST,"ls $zone2011 |");
@dirs=<PLIST>;
close(PLIST);

open(CONDOR,">for_condor.txt");
print(CONDOR "Executable = ${exe}\n");
print(CONDOR "Universe = vanilla\n");
print(CONDOR "Requirements = Memory > 400  && (Arch==\"X86_64\") && (SlotId!=5 && SlotId!=11)\n");
print(CONDOR  "+CondorGroup=\"cmsfarm\"\n");

foreach $file11 (@dirs) {
    chomp $file11;
    next if (!($file11=~/MW([0-9]+)-MN([0-9]+)/));
    $mw=$1; $mn=$2;
    $file11=$zone2011."/".$file11;

    $file10=$zone2010."/MW${mw}-MN${mn}.root";
    if (! (-e $file10)) {
	print "2011 only $mw $mn\n";
	
	$ofname="$fileLoc/limit_".$mw."_".$mn;
	$cmd="./makeLimitFileELEC2011.exe $lumi2011 $mw $file11 $file11 $ofname";
    } else {
	print "2010/2011 $mw $mn\n";

	$ofname="$fileLoc/limit_".$mw."_".$mn;
	$cmd="./makeLimitFileELEC2Y.exe $mw $lumi2011 $file11 $file11 $file10 $file10 $ofname";
    }
    
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys}\n";
    print CONDOR "Queue \n";
    
}
    
close(CONDOR);
system("condor_submit for_condor.txt");


