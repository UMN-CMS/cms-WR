#!/usr/bin/perl


$a=<>;

@sections=split(/:|\]./,$a);

$run=0;
foreach (@sections) {
    if (/\"([0-9]+)\"/) {
	$run=$1;
    } 
    elsif (/([0-9]+), ([0-9]+)/) {
	print("   '$run:$1-$run:$2',\n");
    } else  { print "$line\n"; }
}


