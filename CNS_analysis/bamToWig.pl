#!perl -w
use strict;

print "track type=wiggle_0 name=coverage description=coverage\n";

my $lastC = "";
while (<>) {
#	print "$_";
 	my @spl = split("\\s+", $_);
	if ($spl[0] ne $lastC) {
		print "variableStep chrom=$spl[0]\n";
	}
 	$lastC=$spl[0];
 	my $depth = 0;
 	for( my $i=3; $i< scalar(@spl) ; $i=$i+3  ){
 		$depth = $depth + $spl[$i];
 	}
	print "$spl[1]\t$depth\n";
}
