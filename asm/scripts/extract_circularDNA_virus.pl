#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Spec;

my $sfetch = `which esl-sfetch`;
chomp($sfetch);
unless ( $sfetch && -x $sfetch ) {
    die("make sure you have loaded hmmer/3 so that esl-sfetch in path");
}

my $dir = 'blobOut';
my $query_col = 5;
my $query = 'Viruses';
my $genomedir = 'genomes';
my $genomeext = 'vecscreen_shovill.fasta';

my $outdir = 'virus_search';
mkdir($outdir) unless -d $outdir;
my $outseq = Bio::SeqIO->new(-format => 'fasta',
			   -file   => ">".File::Spec->catfile($outdir,
							      'virus_DNA_hits.fasta'));
opendir(my $od => $dir) || die $!;
foreach my $file ( readdir($od) ) {
    next unless ($file =~ /(\S+)\.AA\.blobDB\.table\.txt/);
    my $stem = $1;
    open(my $fh => File::Spec->catfile($dir,$file)) || die $!;
    while(<$fh>) {
	next if /^\#/;
	chomp;
	my @row = split(/\t/,$_);
	my ($ctgname,$target) = ($row[0],$row[$query_col]);
	if ( $target =~ /$query/) {
	    print join("\t",$stem,$ctgname,$target),"\n";
	    my $genomefile = File::Spec->catfile($genomedir,
						 sprintf("%s.%s",$stem,$genomeext));
	    if ( ! -f "$genomefile.ssi") {
		`$sfetch --index '$genomefile'`;
	    }
	    open(my $fasta => "$sfetch '$genomefile' $ctgname |") || die $!;
	    my $seqin = Bio::SeqIO->new(-fh => $fasta,
					-format => 'fasta');
	    while ( my $seq = $seqin->next_seq ) {
		my $st = $stem;
		$st =~ s/[\(\)]/_/g;
		$seq->display_id(sprintf("%s_%s",$st,$seq->display_id));
		$outseq->write_seq($seq);
	    }		
	}
    }
}

