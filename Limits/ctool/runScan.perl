#!/usr/bin/perl

use Getopt::Long;

$toys=50;
$step=100;
$method="HybridNew";
$jobBase="default";
$special="";
$systmode=0;
$xsec=0.01;
$minmwAllowed=700;
$maxmwAllowed=3000;
$max_mllqq=-1;
$year=2012;
$channel="mu";

GetOptions("toys=i" => \$toys,
	   "systmode=i" => \$systmode,
	   "method=s" => \$method,
	   "minmw=i" => \$minmwAllowed,
	   "maxmw=i" => \$maxmwAllowed,
	   "xsec=f" => \$xsec,
	   "channel=s" => \$channel,
	   "step=i" => \$step,
	   "year=i" => \$year,
	   "maxllqq=i" => \$max_mllqq,
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
$lumi2011=4990;


#$lumi2011=4700;
if ($channel eq "e-mu") {
    $ratesdb="ratesdb.csv,ratesdb_elec.csv";
    
    $data2012="/local/cms/user/dahmes/wr2012/HPAResults/GoodRuns/run2012AB/data-mu-top-qcd-ichep-jun22.root".",".
	"/local/cms/user/franzoni/WR2012/JUN22bis.analysis.V02-02-00.june7.455files/JUN22bis.analysis.V02-02-00.june7.all_HADDED.root";
    $systdb="systematicsdb_mu_2012.csv,systematicsdb_elec_2012.csv";
    $lumi2012="3600,3600";
    $xsec=$xsec.",".$xsec;
} elsif ($channel eq "mu2ecm") {
    $ratesdb="ratesdb.csv,ratesdb.csv";
    
    $data_both=$data2011.",".
	"/local/cms/user/dahmes/wr2012/HPAResults/GoodRuns/run2012AB/data-mu-top-qcd-ichep-jun22.root";
    $systdb="systematicsdb.csv,systematicsdb_elec_2012.csv";
    $lumi2012=3600;
    $xsec=$xsec;
} elsif ($channel=~/mu/) {
    $ratesdb="ratesdb.csv";
    $data2012="/local/cms/user/dahmes/wr2012/HPAResults/GoodRuns/run2012AB/data-mu-top-qcd-ichep-jun22.root";
    $lumi2012=3600;

    if ($year==2012) {
	$systdb="systematicsdb_mu_2012.csv";
    } else {
	$systdb="systematicsdb.csv";
    }

} else {
    $ratesdb="ratesdb_elec.csv";
    $data2012="/local/cms/user/franzoni/WR2012/JUN22bis.analysis.V02-02-00.june7.455files/JUN22bis.analysis.V02-02-00.june7.all_HADDED.root";
    $lumi2012=3600;

    if ($year==2012) {
	$systdb="systematicsdb_elec_2012.csv";
    }

}

$fileLoc=$workloc."/input";
system("mkdir -p $fileLoc");
system("rm -rf ${workloc}/roostats*");

$ratesdblist=$ratesdb;
$ratesdblist=~s/,.*//;

open(PLIST,$ratesdblist);
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
    if ($year==2012) {
	next if (!($item=~/_2012/));
    } else {   
	next if (!($item=~/_2011/));
    }

    $mw=$1; $mn=$2;

    next if (($mw%$step)!=0);

    next if ($mw < $minmwAllowed);
    next if ($mw > $maxmwAllowed);

    print "Got $mw $mn ($channel)\n";

    $mwr_mn{$mw}=$mn;

    $mwmin=$mw if ($mwmin>$mw);
    $mwmax=$mw if ($mwwax<$mw);

    $ofname="$fileLoc/limit_".$mw."_".$mn;
    if ($channel eq "mu2ecm") {
	$mn2011=(int($mw/200)+1)*100;
	$mntext=$mn2011.",".$mn;
	$xseceff=$xsec*xsecratio($mw);
	    	
	$mode=$channel;
	$cmd="./makeLimitFile.exe -m $mode -l $lumi2011,$lumi2012 -w $mw -n $mntext -y 2011,2012 -x $xseceff,$xsec -d $data_both -r $ratesdb -o $ofname -s $systdb ";

	$comments="xsec=$xsec,$xseceff";
    } elsif ($year==2011) {
	$cmd="./makeLimitFile.exe -l $lumi2011 -w $mw -n $mn -x $xsec -d $data2011 -r $ratesdb -o $ofname -s $systdb ";
	$comments="xsec=$xsec";

    } else {

	$mode=$channel;
	$cmd="./makeLimitFile.exe -m $mode -l $lumi2012 -w $mw -n $mn -y $year -x $xsec -d $data2012 -r $ratesdb -o $ofname -s $systdb ";
	if ($max_mllqq>0) {
	    $cmd.=" -b $max_mllqq";
	}

	$comments="xsec=$xsec";
    }
    system($cmd);

    $mass=sprintf("%04d%04d",$mw,$mn);


    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${mw}_${mn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
    print CONDOR "Queue \n";
    
}
# interpolations
#for ($amw=$mwmin; $amw<$mwmax; $amw+=$step) {
#    next if (($amw%100)==0);
#    $mwb=((int($amw/100))*100);
#    $mwa=((int($amw/100)+1)*100);
#    $mnb=$mwr_mn{$mwb};
#    $mna=$mwr_mn{$mwa};
##    print "$mwa $mwb $mna $mnb\n";
#
#    $amn=$mnb;
#	
#    $ofname="$fileLoc/limit_".$amw."_".$amn;
#    $cmd="./makeLimitFile.exe -l $lumi2011 -w $amw -n $amn -x $xsec -d $data2011 -o $ofname -r $ratesdb -s $systdb ";
#    $cmd.=sprintf(" -I %d,%d:%d,%d",$mwb,$mnb,$mwa,$mna);
#    print "$cmd\n";
#    system($cmd);
#
#    $mass=sprintf("%04d%04d",$amw,$amn);
#
#    $comments="xsec=$xsec";
#
#    print CONDOR "Arguments = ${pwd} ${workloc} ${ofname} ${mass} ${workloc}/limit_${amw}_${amn}.log ${method} ${toys} \\\"${comments}\\\" ${special}\n";
#    print CONDOR "Queue \n";
#    
#}
 

   
close(CONDOR);
#system("condor_submit for_condor.txt");


sub xsecratio($) {

    @ratios=(1.346322,
	     1.380845,
	     1.401951,
	     1.447967,
	     1.471503,
	     1.516183,
	     1.574026,
	     1.632000,
	     1.683206,
	     1.724256,
	     1.829525,
	     1.862078,
	     1.956977,
	     2.068023,
	     2.137391,
	     2.221935,
	     2.361538,
	     2.478674,
	     2.605128,
	     2.768387,
	     2.923077,
	     3.025788,
	     3.269565,
	     3.455446);
    
    my($amw)=@_;
    return 1.0/$ratios[($amw/100)-7];    
}
