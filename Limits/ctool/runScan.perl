#!/usr/bin/perl

use Getopt::Long;

$toys=50;
$step=100;
$method="HybridNew";
$jobBase="default";
$special="";
$systmode=0;
$xsec=0.01;
$minmwAllowed=1000;
$maxmwAllowed=2500;

GetOptions("toys=i" => \$toys,
	   "systmode=i" => \$systmode,
	   "method=s" => \$method,
	   "minmw=i" => \$minmwAllowed,
	   "maxmw=i" => \$maxmwAllowed,
	   "xsec=f" => \$xsec,
	   "step=i" => \$step,
	   "special=s" => \$special,
	   "jobname=s" => \$jobBase);



$pwd=`pwd`;
chomp($pwd);
$exe=$pwd."/scanPoint.sh";
$workloc="/local/cms/user/".$ENV{"USER"}."/heavyNuLimits";
if ($jobBase ne "default") {
    $workloc=$workloc."/${jobBase}";
}

$data2011="/local/cms/user/dahmes/wr2011/data_run2011A_run2011B/data-run2011a-run2011b-dec23.root";

$lumi2011=4962;
$systdb="systematicsdb.csv";
$ratesdb="ratesdb.csv";

$fileLoc=$workloc."/input";
system("mkdir -p $fileLoc");
system("rm -rf ${workloc}/roostats*");

open(PLIST,$ratesdb);
@items=<PLIST>;
close(PLIST);

open(CONDOR,">for_condor.txt");
print(CONDOR "Executable = ${exe}\n");
print(CONDOR "Universe = vanilla\n");
print(CONDOR "Requirements = Memory > 400  && (Arch==\"X86_64\") ");

#print(CONDOR "&& (SlotId>=5 && SlotId<=10)");
#print(CONDOR "&& (SlotId==6 || SlotId==11)");
print(CONDOR "&& (SlotId==6)");
print(CONDOR "\n");
print(CONDOR  "+CondorGroup=\"cmsfarm\"\n");

$mwmin=10000;
$mwmax=0;
%mwr_mn={};

# direct mass points 
foreach $item (@items) {
    chomp $item;
#    print $item;
#    next if (!($item=~/WRToNuLeptonToLLJJ_MW-([0-9]+)_MNu-([0-9]+)/));
    next if (!($item=~/signal_([0-9]+)_([0-9]+)/));
    $mw=$1; $mn=$2;

    next if ($mw < $minmwAllowed);
    next if ($mw > $maxmwAllowed);

    $mwr_mn{$mw}=$mn;

    $mwmin=$mw if ($mwmin>$mw);
    $mwmax=$mw if ($mwwax<$mw);

    $ofname="$fileLoc/limit_".$mw."_".$mn;
    $cmd="./makeLimitFile.exe -l $lumi2011 -w $mw -n $mn -x $xsec -d $data2011 -r $ratesdb -o $ofname -s $systdb ";
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);

    $comments="xsec=$xsec";

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";
    
}
# interpolations
for ($amw=$mwmin; $amw<$mwmax; $amw+=$step) {
    next if (($amw%100)==0);
    $mwb=((int($amw/100))*100);
    $mwa=((int($amw/100)+1)*100);
    $mnb=$mwr_mn{$mwb};
    $mna=$mwr_mn{$mwa};
#    print "$mwa $mwb $mna $mnb\n";

    $amn=$mnb;
	
    $ofname="$fileLoc/limit_".$amw."_".$amn;
    $cmd="./makeLimitFile.exe -l $lumi2011 -w $amw -n $amn -x $xsec -d $data2011 -o $ofname -r $ratesdb -s $systdb ";
    $cmd.=sprintf(" -I %d,%d:%d,%d",$mwb,$mnb,$mwa,$mna);
    print "$cmd\n";
    system($cmd);

    $mass=sprintf("%04d%04d",$amw,$amn);

    $comments="xsec=$xsec";

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${amw}_${amn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";
    
}
 

   
close(CONDOR);
system("condor_submit for_condor.txt");

