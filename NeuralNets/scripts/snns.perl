#!/usr/bin/perl

#
# This application runs batchman for training SNNS
# neural networks.
#



# Configuration
###############
$checkout="nnbackeff";
$checkoutllr="nnllr";
$snnsTools="/home/jmmans/src/SNNS/SNNSv4.3/tools/bin/i686-pc-linux-gnu";

# defaults
##

$sourcenetdef          =   "blank_network.net";
$cutoffdef             =   "0.9";
$cyclesdef             =   "400";

@backgrounds = ("madgraph_zmumujets_m180_reco385_nntraining.txt",
		"pythia6_tauola_ttbar_mumu_tuneZ2_reco385_nntraining.txt");


###############
# End Config

@INC=(".",@INC);

require jobstools;

print <<MENU;

SNNS BatchMan Script Tool

MENU

    @masspoints=();
open(SRCS,"ls ../patterns/heavyNuReco_WR*|");

while (<SRCS>) {
    if (/WR([0-9]+)_nuRmu([0-9]+)/) {
	push @masspoints,($1*10000)+$2;
    }
}
close(SRCS);

foreach $masspoint (@masspoints) {
    $mw=sprintf("%d",$masspoint/10000);
    $mnu=$masspoint%10000;
    print "Working on $mw, $mnu ...\n";

    makePatternFiles($mw,$mnu);
    doMassPoint($mw,$mnu);

    system("$snnsTools/batchman -f /tmp/snns.bat -q");

    system("$snnsTools/snns2c $outNet nn_${mw}_${mnu}.c heavyNuNN_${mw}_${mnu}");
    system("mv nn_${mw}_${mnu}.c ../src/");
    unlink("nn_${mw}_${mnu}.h");

    unlink($trainingFileCreated);
    unlink($testingFileCreated);
}

open(SRC,">../src/HeavyNuNetworks.cc");
print SRC "#include \"HeavyNu/NeuralNets/interface/HeavyNuNetworks.h\"\n\nextern \"C\" {\n";

foreach $masspoint (@masspoints) {
    $mw=sprintf("%d",$masspoint/10000);
    $mnu=$masspoint%10000;
    print SRC "  int heavyNuNN_${mw}_${mnu}(float* in, float* out, int init);\n";
}
   
print SRC <<"EOG";
}

float HeavyNuNetworks::evaluate(int mw, int mnu, const std::vector<float>& invalstd) {
    float retval=-1000f;
    float inval[2]; 
    inval[0]=invalstd[0];    inval[1]=invalstd[1];
    
EOG

$i=0;
foreach $masspoint (@masspoints) {
    $mw=sprintf("%d",$masspoint/10000);
    $mnu=$masspoint%10000;
    print SRC "  ";
    print SRC "else " if ($i!=0);
    print SRC "if (mw==${mw} && mnu==${mnu}) heavyNuNN_${mw}_${mnu}(inval, &retval, 0);\n";
    $i++;
}

print SRC "  return retval;\n}\n";
close(SRC);


exit;

sub makePatternFiles {
    my($mw,$mnu)=@_;

    my @trainLines=();
    my @testLines=();

    $signalPats=sprintf("../patterns/heavyNuReco_WR%d_nuRmu%d_1-anal_nntraining.txt",$mw,$mnu);
    open(SIGNAL,$signalPats) || die "No such file $signalPats";
    while (<SIGNAL>) {
	$l1=$_;
	$l2=<SIGNAL>;
	if (int(rand(2))==1) {
	    push @trainLines,$l1,$l2;
	} else {
	    push @testLines,$l1,$l2;
	}
    }
    close(SIGNAL);

    foreach $backFile (@backgrounds) {
	open(BACK,"../patterns/$backFile");
	while (<BACK>) {
	    $l1=$_;
	    $l2=<BACK>;
	    if (int(rand(2))==1) {
		push @trainLines,$l1,$l2;
	    } else {
		push @testLines,$l1,$l2;
	    }
	}
	close(BACK);	
    }

    my $ntrain=int($#trainLines/2+0.5);
    my $ntest=int($#testLines/2+0.5);

    $trainingFileCreated="/tmp/training_".$mw."_".$mnu.".pat";
    open(TRFILE,">$trainingFileCreated");
    print(TRFILE "SNNS pattern definition file V3.2\n");
    print(TRFILE "generated at Thu Apr 21 19:56:37 1994\n\n\n");
    print(TRFILE "No. of patterns : $ntrain\n");
    print(TRFILE "No. of input units : 2\n");
    print(TRFILE "No. of output units : 1\n\n");

    foreach $line (@trainLines) {
	print(TRFILE $line);
    }
    close(TRFILE);

    $testingFileCreated="/tmp/testing_".$mw."_".$mnu.".pat";
    open(TSTFILE,">$testingFileCreated");
    print(TSTFILE "SNNS pattern definition file V3.2\n");
    print(TSTFILE "generated at Thu Apr 21 19:56:37 1994\n\n\n");
    print(TSTFILE "No. of patterns : $ntest\n");
    print(TSTFILE "No. of input units : 2\n");
    print(TSTFILE "No. of output units : 1\n\n");

    foreach $line (@testLines) {
	print(TSTFILE $line);
    }
    close(TSTFILE);
    
#    open(TSTFILE,">/tmp/testing_$mw_$mnu.pat");

    

    
}

sub doMassPoint {
    my($mw,$mnu)=@_;
    my($sourcenet, $trainingpat, $validationpat, $cutoff, $resfile);
    $resfile="/tmp/snns.res";
    $useLLR="n";
#    if (jobstools::prompt("  Change options?","n")=~/[yY]/) {
#	$sourcenet=jobstools::prompt("  Source network?",$sourcenetdef{$nc});
#	$trainingpat=jobstools::prompt("  Training pattern?",$trainingpatdef{$nc});
#	$validationpat=jobstools::prompt("  Testing pattern?",$validationpatdef{$nc});
#	$outNet=$sourcenet;
#	$cycles=jobstools::prompt("  Cycles?",$cyclesdef{$nc});
#	$useLLR=jobstools::prompt("  Use LLR?","n");
#	if ($useLLR=~/[yY]/) {
#	    $sw=jobstools::prompt("    Signal weight?","1.0");
#	    $bw=jobstools::prompt("    Background weight?","1.0");
#	} else {
#	    $cutoff=jobstools::prompt("  Efficiency cutoff?",$cutoffdef{$nc});
#	}
#    } else {
    $sourcenet=$sourcenetdef;
    $trainingpat=$trainingFileCreated;
    $validationpat=$testingFileCreated;
    $outNet="snns_".$mw."_".$mnu.".net";
    $cycles=$cyclesdef;
    $cutoff=$cutoffdef;
#    }    
    open(BATCHFILE,">/tmp/snns.bat") || die "Could not open!";
    if ($useLLR=~/[yY]/) {
	makebatchllr($sourcenet, $trainingpat, $validationpat, $outNet, $sw, $bw, $resfile,$cycles);
    } else {
	makebatchman($sourcenet, $trainingpat, $validationpat, $outNet, $cutoff, $resfile,$cycles);
    }
    close(BATCHFILE);
}


sub makebatchman {
    my($sourcenet, $trainingpat, $validationpat, $outNet, $cutoff, $resfile, $maxcycles)=@_;
       print BATCHFILE <<EOBATCH;
testcmd:="$checkout $resfile $cutoff";
oldnbkg:= 100000;
nbkg:= 100000;

loadNet("$sourcenet");
loadPattern("$trainingpat");
loadPattern("$validationpat");
setInitFunc("Randomize_Weights",1.0,-1.0);
setLearnFunc("Rprop",0.2,50,5);
initNet();

setPattern("$trainingpat");
#setClassDistrib(TRUE,1,20)

setSeed();
setShuffle(TRUE);

while CYCLES < $maxcycles and SIGNAL == 0 do
  trainNet()
  if CYCLES mod 2 == 0 then
    setPattern("$validationpat");
    testNet()
    saveResult("$resfile",1,PAT,FALSE,TRUE,"create")
    execute(testcmd,nbkg)
    if (nbkg<oldnbkg) then
      saveNet("$outNet");
      oldnbkg:=nbkg;
      print("(",CYCLES,") Saving at ",nbkg," background events...")
    else
      print("(",CYCLES,") Not saving at ",nbkg," background events...")
    endif
    setPattern("$trainingpat");
  endif
endwhile

execute("rm $resfile");
EOBATCH
}

sub makebatchllr {
    my($sourcenet, $trainingpat, $validationpat, $outNet, $sigw, $backw, $resfile, $maxcycles)=@_;
       print BATCHFILE <<EOBATCH2;
testcmd:="$checkoutllr $resfile $sigw $backw";
oldllr:= -1;
llr:= 100000;

loadNet("$sourcenet");
loadPattern("$trainingpat");
loadPattern("$validationpat");
setInitFunc("Randomize_Weights",1.0,-1.0);
setLearnFunc("Rprop",0.2,50,5);
initNet();

setPattern("$trainingpat");
#setClassDistrib(TRUE,2,1)

setSeed();
setShuffle(TRUE);

while CYCLES < $maxcycles and SIGNAL == 0 do
  trainNet()
  if CYCLES mod 2 == 0 then
    setPattern("$validationpat");
    testNet()
    saveResult("$resfile",1,PAT,FALSE,TRUE,"create")
    execute(testcmd,llr)
    if (llr>oldllr) then
      saveNet("$outNet");
      oldllr:=llr;
      print("(",CYCLES,") Saving at llr=",llr)
    else
      print("(",CYCLES,") Not saving at llr=",llr)
    endif
    setPattern("$trainingpat");
  endif
endwhile

execute("rm $resfile");
EOBATCH2
}



sub selectlist {
    local @thelist=sort(@_);
    local $option=1;

    do {
	$option=1;
	foreach (@thelist) {
	    if (/[^.]\w/) {
		print "$option ) $_\n";
		$option ++;
	    }
	}
	print "\n   Select : ";
	$option=<STDIN>;
	$option=~s/[ ]+/,/g;
	$option=~s/\s//g;
    } while ($option=~/[^0-9 ,]/);
    @alist=split(/[ ,]+/,$option);
    @retlist=();
    foreach (@alist) {
	push(@retlist,$thelist[$_+1]);
#	print $thelist[$_+1];
    }

    return @retlist;
}
