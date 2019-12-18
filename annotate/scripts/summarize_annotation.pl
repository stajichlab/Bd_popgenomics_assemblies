#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;

my $strains = 'samples.tab';
my $dir = 'funannotate';

my %strain_stats;
my %stats;
my %header_seen;

open(my $fh => $strains) || die $!;
while(<$fh>) {
    my ($strain) = split;
    $strain =~ s/[\(\)]//;

    my $predfolder = "$dir/$strain/annotate_results";
    if ( -d $predfolder ) {
	opendir(my $folder => $predfolder) || die "cannot read $predfolder: $!";
	foreach my $resultfile ( readdir($folder) ) {
	    my $fullfile = File::Spec->catfile($predfolder,$resultfile);
	    if ( $resultfile =~ /\.gff3$/ ) {
		my $gff = $fullfile;
		open(my $gffin => $gff ) || die $!;
		while(<$gffin> ) {
		    chomp;
		    last if /^\#{1,2}\s+FASTA/;
		    next if /^\#/;
		    my @row = split(/\t/,$_);
		    my %grp = (map { split(/=/,$_) } split(/;/,pop @row));
		    my $group = exists $grp{Parent} ? $grp{Parent} : $grp{ID};
		    $header_seen{$row[2]} = 1;
		    $stats{$strain}->{$row[2]}->{$group}++
		}
	    }
	}
    }
}

#my @header = sort keys %header_seen;
# force order
my @header;
# here we force the expected order of features, skip a feature if it wasn't seen in the
# dataset
for my $h (qw(gene mRNA CDS exon five_prime_UTR three_prime_UTR tRNA) ) {
    if ( exists $header_seen{$h} ) {
	push @header, $h;
	delete $header_seen{$h};
    }
}
# and add in any other features types that are seen but weren't in the list above
push @header, sort keys %header_seen;

for my $strain ( keys %stats) {
    for my $type ( @header ) {
	$strain_stats{$strain}->{$type} = scalar keys %{$stats{$strain}->{$type}};
    }
}
print join("\t", qw(SampleID), @header), "\n";
foreach my $sp ( sort keys %stats ) {
    print join("\t", $sp, map { $strain_stats{$sp}->{$_} || 'NA' } @header), "\n";
}
