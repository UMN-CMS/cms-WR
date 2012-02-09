#!/usr/bin/perl

use Getopt::Long;

$toys=50;
$method="HybridNew";
$jobBase="default";
$special="";
$systmode=0;
$xsec=0.01;

GetOptions("toys=i" => \$toys,
	   "systmode=i" => \$systmode,
	   "method=s" => \$method,
	   "xsec=f" => \$xsec,
	   "special=s" => \$special,
	   "jobname=s" => \$jobBase);



$pwd=`pwd`;
chomp($pwd);
$exe=$pwd."/scanPoint.sh";
$workloc="/local/cms/user/".$ENV{"USER"}."/heavyNuLimits";
if ($jobBase ne "default") {
    $workloc=$workloc."/${jobBase}";
}

#$data2010="/local/cms/user/jmmans/heavyNu_jun14/mu/run2010_42X_june4.root";
#$signal2010="/local/cms/user/dudero/heavyNuAnalysis/dirTrackRelIso+MinDR=0.5/dirSignal";
#$signal2010="/local/cms/user/jmmans/heavyNu_jun14/mu";

#$data2011="/local/cms/user/dahmes/wr2011/MuResults/GoodRuns/may13-191ipb/heavyNuAnalysis/data2011_191ipb.root";
#$data2011="/local/cms/user/jmmans/heavyNu_oct2x/mu/run2011A-nominal_2375_oct25.root";
$data2011="/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root";
$signal2011="/local/cms/user/jmmans/heavyNu_signal/";


$lumi2011=4700;
$systdb="systematicsdb.csv";

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
print(CONDOR "Requirements = Memory > 400  && (Arch==\"X86_64\") ");

#print(CONDOR "&& (SlotId>=5 && SlotId<=10)");
print(CONDOR "&& (SlotId==6 || SlotId==11)");
print(CONDOR "\n");
print(CONDOR  "+CondorGroup=\"cmsfarm\"\n");
 
foreach $item (@items) {
    chomp $item;
#    print $item;
#    next if (!($item=~/WRToNuLeptonToLLJJ_MW-([0-9]+)_MNu-([0-9]+)/));
    next if (!($item=~/signal_([0-9]+)_([0-9]+).root/));
    $mw=$1; $mn=$2;

    $s11=$signal2011."/".$item;

    if (! -e $s11) {
	print "No $s11\n";
	next;
    }

	
    $ofname="$fileLoc/limit_".$mw."_".$mn;
    $cmd="./makeLimitFile.exe -l $lumi2011 -w $mw -n $mn -x $xsec -d $data2011 -i $s11 -o $ofname -s $systdb ";
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);

    $comments="xsec=$xsec";

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";
    
}
    
close(CONDOR);
system("condor_submit for_condor.txt");

