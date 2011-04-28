#!/usr/bin/perl
open (IN, 'Shermin/shermin.txt') || die ("can't open shermin.txt\n");
@db = <IN>;
close IN;

#opendir(IMD, "reference_sequences") || die("Cannot open directory");
#@dir = readdir(IMD);
#closedir IMD;

@interest;
foreach $l (@db) {
	@line = split(/\|/, $l);
#	foreach (@line) {
#		print $_ . "\n";
#	}
	
	
#	foreach $f (@dir) {
		$fna = "";
		@line[0] =~ /(^NC_\w+)\s*/;
		$file = "reference_sequences/$1.fna.fna";
		#print "$file \n";
		if ( -e $file) {
		#print ("$f @line[0]\n");
			open(FILE, "$file") || die ("can't open NC file");
			while ( $curline = <FILE>) {
				if ( $curline !~ /\>/) {
					chomp $curline;
					$fna .= $curline;
				}
				if (length($fna) >= @line[3]) {
					last;
				}
			}
			$seqInterest = substr($fna, @line[3]-150, 150);
			$seqInterest = ">$file\n" . $seqInterest;
			push(@interest, $seqInterest );
		}
#	}
}
foreach $l (@interest) {
	print "$l\n";
}