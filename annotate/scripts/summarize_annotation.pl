#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;

my $strains = 'samples.dat';
my $dir = 'funannotate';

my %stats;
my @header;
my %header_seen;

open(my $fh => $strains) || die $!;
while(<$fh>) {
    my ($srr,$strain) = split;
    $strain =~ s/[\(\)]//;
    
    # protein count
    my $name = 'Proteins';
    $stats{$strain}->{$name} = 0;
    if( ! exists $header_seen{$name} ) {
	push @header, $name;
	$header_seen{$name} = 1;
    }
    my $predfolder = "$dir/$strain/predict_results";
    if ( -d $predfolder ) {
	opendir(my $folder => $predfolder) || die "cannot read $predfolder: $!";
	foreach my $resultfile ( readdir($folder) ) {
	    next unless $resultfile =~ /\.proteins.fa$/;
	    my $pepfile = File::Spec->catfile($predfolder,$resultfile);
	    open(my $grp => "grep -c '^>' $pepfile |") || die $!;
	    while(<$grp>) {
		s/^\s+//;
		my ($n) = split;
		$stats{$strain}->{$name} = $n;
	    }
	    
	}
    }
}

print join("\t", qw(SampleID), @header), "\n";
foreach my $sp ( sort keys %stats ) {
    print join("\t", $sp, map { $stats{$sp}->{$_} || 'NA' } @header), "\n";
}
