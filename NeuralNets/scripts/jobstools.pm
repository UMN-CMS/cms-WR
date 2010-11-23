package jobstools;
#
#
# general tools....
#
#


sub printend {
    @endinglist=(@endinglist,@_);
}

sub printl {
    @outputlist=(@outputlist,@_);
}

sub writeparams {
    my($output)=@_;
    open(PARAMETER,">$output");
    foreach(@outputlist) {
	print(PARAMETER $_);
    }
    foreach(@endinglist) {
	print(PARAMETER $_);
    }
    close(PARAMETER);
}

sub prompt($$) {
    local ($prompt,$default)=@_;
    local $local;
    print("$prompt [$default] : ");
    $local=<STDIN>;
    chop $local;
    if (length($local)==0) { $local=$default; }
    return $local;
}

1;
