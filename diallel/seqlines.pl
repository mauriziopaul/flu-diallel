#!/usr/bin/perl
use English;
use FileHandle;
use integer;

unless (@ARGV){
    print STDERR "seqlines from to by infile\n"; exit 1;
}

my($from, $to, $by, $inFile) = @ARGV;

my $ist = FileHandle->new($inFile, "r");
if (not defined $ist){
	die "Could not open $inFile";
}
my $i = 0;
while (my $line = $ist->getline){
    $i++;
    if ($to ne "END"){
    	if ($i > $to){
    		last;
    	}
    }
    if ($i < $from) {
    	next;
    }
    #print( ($i - $from + 1) % $by, "\n");
    if (0 == (($i - $from) % $by)) {
        print STDOUT $line;
    }
}
$ist->close;
