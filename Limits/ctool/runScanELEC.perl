#!/usr/bin/perl

use Getopt::Long;

$toys=50;
$interpol=100;
$method="HybridNew";
$jobBase="default";
$special="";
$systmode=0;
$xsec=0.01;

GetOptions("toys=i" => \$toys,
	   "systmode=i" => \$systmode,
	   "method=s" => \$method,
	   "xsec=f" => \$xsec,
	   "interpol=i" => \$interpol,
	   "special=s" => \$special,
	   "jobname=s" => \$jobBase);



$pwd=`pwd`;
chomp($pwd);
$exe=$pwd."/scanPoint.sh";
$workloc="/local/cms/user/".$ENV{"USER"}."/heavyNuLimits";
if ($jobBase ne "default") {
    $workloc=$workloc."/${jobBase}";
}

$signal2011="/local/cms/user/jmmans/heavyNuINR/2011_final/";

$lumi2011=4700;
$systdb="systematicsdbELEC.csv";

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

$mwmin=10000;
$mwmax=0;
%mwr_mn={};

# direct mass points 
foreach $item (@items) {
    chomp $item;
#    print $item;
#    next if (!($item=~/WRToNuLeptonToLLJJ_MW-([0-9]+)_MNu-([0-9]+)/));
    next if (!($item=~/MW([0-9]+)-MN([0-9]+).root/));
    $mw=$1; $mn=$2;

    $mwr_mn{$mw}=$mn;

    $mwmin=$mw if ($mwmin>$mw);
    $mwmax=$mw if ($mwwax<$mw);

    $s11=$signal2011."/".$item;

    if (! -e $s11) {
	print "No $s11\n";
	next;
    } 
     

	
    $ofname="$fileLoc/limit_".$mw."_".$mn;
    $cmd="./makeLimitFileELEC.exe -l $lumi2011 -w $mw -n $mn -x $xsec -d $s11 -i $s11 -o $ofname -s $systdb ";
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);

    $comments="xsec=$xsec";

    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";
    
}
# interpolations
for ($amw=$mwmin; $amw<$mwmax; $amw+=$interpol) {
    next if (($amw%100)==0);
    $mwb=((int($amw/100))*100);
    $mwa=((int($amw/100)+1)*100);
    $mnb=$mwr_mn{$mwb};
    $mna=$mwr_mn{$mwa};
    print "$mwa $mwb $mna $mnb\n";


    $fa=sprintf("%s/MW%d-MN%d.root",$signal2011,$mwa,$mna);
    $fb=sprintf("%s/MW%d-MN%d.root",$signal2011,$mwb,$mnb);

    if (! -e $fa) {
	print "No $fa\n";
	next;
    }
    if (! -e $fb) {
	print "No $fb\n";
	next;
    }

    $amn=$mnb;
	
    $ofname="$fileLoc/limit_".$amw."_".$amn;
    $cmd="./makeLimitFile.exe -l $lumi2011 -w $amw -n $amn -x $xsec -d $s11 -o $ofname -s $systdb ";
    $cmd.=sprintf(" -I %d,%d,%s,%d,%d,%s",$mwb,$mnb,$fb,$mwa,$mna,$fa);
    print "$cmd\n";
    system($cmd);


    $mass=sprintf("%04d%04d",$amw,$amn);
    
    $comments="xsec=$xsec";
    
    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${amw}_${amn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";


}
 

   
close(CONDOR);
#system("condor_submit for_condor.txt");

