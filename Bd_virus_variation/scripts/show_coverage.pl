#!/usr/bin/env perl
use strict;
use warnings;


my $genome_cov = "coverage";

my %cov;
my %chroms;
foreach my $dir ( $genome_cov ) {
    opendir(my $gcov => $dir) || die "cannot open $dir: $!";
    for my $file ( readdir($gcov) ) {
    	next unless ( $file =~ /(\S+)\.regions\.bed\.gz$/);
	my $strain = $1;
	$strain =~ s/[\(\)]//g;
    	open(my $fh => "gzip -dc \"$dir/$file\" |");
	while(<$fh>) {
	    my ($chrom,$start,$end,$coverage) = split;
	    if ( ! exists $chroms{$chrom} ) {
		$chroms{$chrom} = abs($end - $start);
	    }
	    $cov{$strain}->{$chrom} = $coverage;
	}
    }
}

my @chroms = sort keys %chroms;
print join("\t", 'STRAIN',@chroms), "\n";
for my $strain ( sort keys %cov) {
    print join("\t", $strain, map { sprintf("%.2f",$cov{$strain}->{$_} || 0) } @chroms),"\n";
}
