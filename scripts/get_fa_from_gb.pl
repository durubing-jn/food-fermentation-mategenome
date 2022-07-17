#!/usr/bin/perl
use strict;
use warnings;

my $id;
my $seq;

open IN,$ARGV[0] or die "Can't open genebank file!";
while (<IN>) {
	chomp;
	if (/ACCESSION/) {
		my @infor = split (/\s+/,$_);
		$id = $infor[1];
		print ">$id\n";
	}
	if (/\s+\d+\s[ATGCatgc]{10}\s/) {
		$seq = $_;
		$seq =~ s/\s+//g;
		$seq =~ s/\d+//g;
		print "$seq\n";
	}
}
close IN;

