#!/usr/bin/perl

open(MYINPUTFILE, "<foo.txt");
my %lumis  = ();
while(<MYINPUTFILE>) {
    $a=<>;
    @sections=split(':',$a);
    
    $run=0;
    foreach (@sections) {
        if (/\"([0-9]+)/) {
            $run=$1;
            if (exists $lumis{$run}) {
                # print "I already found this run: $run\n";
                # for my $key ( keys %lumis ) {
                #     my @value = @$lumis{$key};
                #     foreach (@{$lumis{$key}}) {
                #         push(@value,$_);
                #     }
                #     print "$key => @value\n";
                # }
            } else {
                # print "ERASING!\n";
                delete $lumis{$run};
                # print "(1) run is $run with size $#{$lumis{$run}}\n";
            }
        } 
        elsif (/([0-9]+)-([0-9]+)\"/) {
            my @values = ();
            # print "(2) run is $run with size $#{$lumis{$run}}\n";
            foreach (@{$lumis{$run}}) {
                push(@values,$_);
            }
            # print "here they are: @values\n";
            for ($ls=$1; $ls<=$2; $ls++) { 
                if ( $run == -1 ) { 
                    print "Trying to add ls $ls to run $run\n";
                }
                # print "$values\n";
#                if (grep {$_ eq $ls} @values) {
#                    print "Lumi section already found: $ls $run\n";
#                } else { 
                    # print "$#values\n";
                    # push(@values,$ls);
                    # print "$#values\n";
                    # print "new: @values\n";
                push(@{$lumis{$run}},$ls);# => \@values;
                if ( $run == -1 ) { 
                    print "Hash content\n";
                    foreach $k (keys %lumis) {
                        print "$k" ; 
                        foreach (@{$lumis{$k}}) {
                            print " : $_";
                        }
                        print "\n";
                    }
                }
                    # for my $key ( keys %lumis ) {
                    #     my @value = $lumis{$key};
                    #     print "$key => @value\n";
                    # }
#                }
            }
            # print("   '$run:$1-$run:$2',\n");
        } # else  { print "$line\n"; }
    }
}

# print "Dump\n";
foreach $run (sort keys %lumis) {
    my @ls = ();
    foreach (@{$lumis{$run}}) {
        push(@ls,$_);
    }
    @ls = sort {$a <=> $b} @ls;
    my $first = shift(@ls);
    my $last = $first;
    my $steps = 0;
    if ( $run == -1 ) { 
        print "Start: $first $last\n" ; 
    }
    while ($#ls >= 0) {
        my $next = shift(@ls);
        if ( $next != $last ) { 
            $steps++;
            if ( $run == -1 ) { 
                print "Interim: $first $last $next $steps $#ls\n";
            }
            if ($next > ($first+$steps)) {
                print "   '$run:$first-$run:$last',\n";
                $first = $next;
                $last = $next;
                $steps = 0;
            } else {
                $last = $next;
            }
        }
    }
    print "   '$run:$first-$run:$last',\n";
}

close(MYINPUTFILE); 

