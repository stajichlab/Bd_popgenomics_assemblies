#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Spec;

my $eslfetch = `which esl-sfetch`;
unless ( -x $eslfetch ) {
    die("make sure you have loaded hmmer/3 so that esl-sfetch in path");
}
my $dir = 'blobOut';
my $query_col = 6;
my $query = 'Viruses';
my $genomedir = 'genomes';
my $genomeext = 'vecclean_shovill.fasta';
my $dir = 'blobOut';
my $outdir = 'virsus_search';
mkdir($outdir) unless -d $outdir;
my $outseq = Bio::Seq->new(-format => 'fasta',
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
	my ($name,$target) = ($row[0],$row[$query_col]);
	print join("\t",$stem,$name,$target),"\n";
	if ( $target =~ /$query/) {
	    my $genomefile = File::Spec->catfile($genomedir,
						 sprintf("%s.%s",$stem,$genomeext));
	    if ( ! -f "$genomefile.ssi") {
		`$eslsfetch --index $genomefile`;
	    }
	    open(my $fasta => "$eslsfetch $genomefile $target |") || die $!;
	    my $seqin = Bio::SeqIO->new(-fh => $fasta,
					-format => 'fasta');
	    while ( my $seq = $seqin->next_seq ) {
		$seq->display_id(sprintf("%s_%s",$stem,$seq->display_id));
		$outseq->write_seq($seq);
	    }		
	}
    }
}

