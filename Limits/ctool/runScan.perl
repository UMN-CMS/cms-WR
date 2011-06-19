#!/usr/bin/perl

use Getopt::Long;

$toys=50;
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

$data2010="/local/cms/user/jmmans/heavyNu_jun14/mu/run2010_42X_june4.root";
#$signal2010="/local/cms/user/dudero/heavyNuAnalysis/dirTrackRelIso+MinDR=0.5/dirSignal";
$signal2010="/local/cms/user/jmmans/heavyNu_jun14/mu";

#$data2011="/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/may13-191ipb/heavyNuAnalysis/data2011_191ipb.root";
$data2011="/local/cms/user/jmmans/heavyNu_jun14/mu/run2011_42X_june4.root";
$signal2011="/local/cms/user/jmmans/heavyNu_jun14/mu";


$lumi2011=204;

$fileLoc=$workloc."/input";
system("mkdir -p $fileLoc");
system("rm -rf ${workloc}/roostats*");

# We use the 2011 as the master list of points
open(PLIST,"ls $signal2011 |");
@items=<PLIST>;
close(PLIST);

open(CONDOR,">for_condor.txt");
print(CONDOR "Executable = ${exe}\n");
print(CONDOR "Universe = vanilla\n");
print(CONDOR "Requirements = Memory > 400  && (Arch==\"X86_64\") && (SlotId>=5 && SlotId<=10)\n");
print(CONDOR  "+CondorGroup=\"cmsfarm\"\n");

foreach $item (@items) {
    chomp $item;
#    print $item;
    next if (!($item=~/2011Sig_mumu_s11_MWR_MNU_([0-9]+)_([0-9]+)/));
    $mw=$1; $mn=$2;

    $s11=$signal2011."/".$item;

    if (! -e $s11) {
	print "No $s11\n";
	next;
    }

    $s10=$signal2010."/heavynu_2010Sig_mumu_s11_MWR_MNU_".$mw."_".$mn.".root";
    if (! (-e $s10)) {
	print "2011 only $mw $mn\n";
	
	$ofname="$fileLoc/limit_".$mw."_".$mn;
	$cmd="./makeLimitFile2011.exe $lumi2011 $mw $data2011 $s11 $ofname";
    } else {
	print "2010/2011 $mw $mn\n";

	$ofname="$fileLoc/limit_".$mw."_".$mn;
	$cmd="./makeLimitFile2Y.exe $mw $lumi2011 $data2011 $s11 $data2010 $s10 $ofname";
    }
    
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys}\n";
    print CONDOR "Queue \n";
    
}
    
close(CONDOR);
#system("condor_submit for_condor.txt");


