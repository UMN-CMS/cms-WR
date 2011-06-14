#!/usr/bin/perl

use Getopt::Long;

$toys=100;
$method="MarkovChainMC";
$jobBase="default";
$mw=1400;
$mn=600;


GetOptions("toys=i" => \$toys,
	   "method=s" => \$method,
	   "jobname=s" => \$jobBase);



$pwd=`pwd`;
chomp($pwd);
$exe=$pwd."/scanSimple.sh";
$workloc="/local/cms/user/".$ENV{"USER"}."/heavyNuSimple";
if ($jobBase ne "default") {
    $workloc=$workloc."/${jobBase}";
}

$data2010="/local/cms/user/dudero/heavyNuAnalysis/dirTrackRelIso+MinDR=0.5/bryansNewMuSkimRun2010A+BcombinedDec22rereco_heavyNuAnalysis.root";
$signal2010="/local/cms/user/dudero/heavyNuAnalysis/dirTrackRelIso+MinDR=0.5/dirSignal";

#$data2011="/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/may13-191ipb/heavyNuAnalysis/data2011_191ipb.root";
$data2011="/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/may17-191ipb-loDiL/heavyNuAnalysis/data2011_loDiL_191ipb.root";
$signal2011="/local/cms/user/pastika/heavyNuAnalysis_2011";

$lumi2011=191;

$fileLoc=$workloc."/input";
system("mkdir -p $fileLoc");
system("rm -rf ${workloc}/roostats*");

# We use the 2011 as the master list of points
open(PLIST,"ls $signal2011 |");
@dirs=<PLIST>;
close(PLIST);

open(CONDOR,">for_condor.txt");
print(CONDOR "Executable = ${exe}\n");
print(CONDOR "Universe = vanilla\n");
print(CONDOR "Requirements = Memory > 400  && (Arch==\"X86_64\") && (SlotId!=5 && SlotId!=11)\n");
print(CONDOR  "+CondorGroup=\"cmsfarm\"\n");

for ($minbin=14; $minbin<50; $minbin++) {

    $s11=$signal2011."/WR".$mw."_nuRmu".$mn."/heavyNuAnalysis_WR".$mw."_nuRmu".$mn.".root";

    if (! -e $s11) {
	print "No $s11\n";
	next;
    }

    $s10=$signal2010."/heavyNuReco_WR".$mw."_nuRmu".$mn."_1-heavyNuAnalysis.root";
    if (! (-e $s10)) {
	$s10=$signal2010."/heavyNuReco_WR".$mw."_nuRmu".$mn."_001-heavyNuAnalysis.root";
    }
    if (! (-e $s10)) {
	next;
	print "2011 only $mw $mn\n";
	
	$ofname="$fileLoc/limit_".$mw."_".$mn;
	$cmd="./makeLimitFile2011.exe $lumi2011 $mw $data2011 $s11 $ofname";
    } else {
	print "2010/2011 $mw $mn\n";

	$ofname="$fileLoc/limit_".$mw."_".$mn."_".$minbin;
	$cmd="./makeLimitFileSimple2Y.exe $mw $minbin $lumi2011 $data2011 $s11 $data2010 $s10 $ofname";
    }
    
    system($cmd);

    $mass=sprintf("%04d%04d%02d",$mw,$mn,$minbin);

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}_${minbin}.log ${method} ${toys}\n";
    print CONDOR "Queue \n";
    
}
    
close(CONDOR);
system("condor_submit for_condor.txt");


