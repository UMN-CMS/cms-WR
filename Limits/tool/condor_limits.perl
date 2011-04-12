#!/usr/bin/perl

#$USER=`whoami`;
#chomp ${USER};
$USER="heavynu";

$dataLoc="/local/cms/user/${USER}/data";
$workLoc="/local/cms/user/${USER}/limits";
$cardLoc="${workLoc}/cards";
$logLoc="${workLoc}/log";
$srcLoc="${workLoc}/src";
$exeLoc="${workLoc}/bin";
$exe="${exeLoc}/batch_doLimitSetting.sh";

$dosyst="yes";

$dataFile="${dataLoc}/data.root";
$lumi=36.1;
$cycles=5000;

system("mkdir -p ${cardLoc}");
system("mkdir -p ${logLoc}");

#Requirements = Memory > 400  && (Arch=="INTEL" || Arch=="X86_64")
open(CONDOR_SCRIPT, ">${workLoc}/condor_script");

print(CONDOR_SCRIPT "Executable = ${exe}\n");
print(CONDOR_SCRIPT "Universe = vanilla\n");
#print(CONDOR_SCRIPT "Requirements = Memory > 400  && (Arch==\"X86_64\") && SlotId!=5 && SlotId!=6 && SlotId!=11 && SlotId!=12 && SlotId!=1 && SlotId!=2\n");
print(CONDOR_SCRIPT "Requirements = Memory > 400  && (Arch==\"X86_64\") && (SlotId!=5 && SlotId!=11)\n");
#print(CONDOR_SCRIPT " && (Machine != \"scorpion22.spa.umn.edu\" && Machine != \"zebra01.spa.umn.edu\")");
#print(CONDOR_SCRIPT " && (Machine != \"zebra02.spa.umn.edu\" && Machine != \"zebra03.spa.umn.edu\")");
print(CONDOR_SCRIPT "+CondorGroup=\"cmsfarm\"\n");

open(FLIST,"ls ${dataLoc}/*.root|");
@flist=<FLIST>;
close(FLIST);

foreach $file (@flist) {
    chomp($file);
    
    $prestub=$file;
    $prestub=~s/_1-anal.root$//;
    $prestub=~s|.*/||;
    next if (!($prestub=~/heavyNu/));

#
# We're setting up a scan over *cross-section values*
#
    $pointsPerDecade=20;
    $xsecpbMin=0.001;
    $xsecpbMax=6.0;

    $decades=log10($xsecpbMax/$xsecpbMin);  # log10(6000) = 3.78 decades
    $totalXsecPts=int(0.5+$pointsPerDecade*$decades); # = 76
    $factor=pow($xsecpbMax/$xsecpbMin,1.0/($totalXsecPts-1)); # = 1.123
    $pointPerJob=2; # minimum allowed by the code, can't go lower!

    $index=0;
    
    for ($pt=0; $pt<$totalXsecPts; $pt+=$pointPerJob) {
	$xl=$xsecpbMin*pow($factor,$pt);
	$xm=$xsecpbMin*pow($factor,$pt+$pointPerJob-1);

	$stub="${prestub}_".sprintf("%02d",$index);
	$cardfile="${workLoc}/cards/${stub}.card";
	$errfile="${workLoc}/log/${stub}.err";
	$outfile="${workLoc}/log/${stub}.log";

	open(CARDFILE,">${cardfile}");
	print CARDFILE "DATA  : ${dataFile}\n";
	print CARDFILE "SIGNALMC : $file\n";
	print CARDFILE "LUMI : ${lumi}\n";
	print CARDFILE "SCAN : LOG $xl $xm $pointPerJob\n";
	print CARDFILE "CYCLES : ${cycles}\n";
	print CARDFILE "DOSYST : ${dosyst}\n";
	print CARDFILE "TUPLEFILE : ${workLoc}/${stub}_limit.root\n";
	close(CARDFILE);

	print CONDOR_SCRIPT "Error = ${errfile}\n";
	print CONDOR_SCRIPT "Output = ${outfile}\n";
	print CONDOR_SCRIPT "Arguments = ${srcLoc} ${cardfile}\n";
	print CONDOR_SCRIPT "Queue \n";
	
	$index++;
    }
}
close(CONDOR_SCRIPT);

#system("condor_submit ${workLoc}/condor_script");


sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

sub pow {
    my $x = shift;
    my $y = shift;
    return exp($y * log($x));
}
