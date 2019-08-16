#!/usr/bin/bash
#SBATCH -p short -n 24 -N 1 --mem 16gb --out logs/viral_search.log
module unload perl
module load parallel
module load fasta
module load prodigal
module load hmmer
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
mkdir -p genomes
pushd genomes
parallel -j $CPU --rpl '{strain}  $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' ln -s {} {strain}.fasta ::: ../asm/*/scaffolds.fasta
popd

OUTBASE=virus_search/gene_search
for QUERY in lib/viral_pep/*.fa
do
    NAME=$(basename $QUERY .fa)
    OUTDIR=$OUTBASE/${NAME}_gene
    mkdir -p $OUTDIR
    parallel -j $CPU  tfastx36 -m 8c -E 1e-5 $QUERY {} \> $OUTDIR/{/.}.TFASTX_rep.tab ::: genomes/*.fasta
done

find $OUTBASE -size 0 | xargs rm

perl scripts/extract_TFASTX_viral_hits.pl

INDIR=virus_search
prodigal -n -i $INDIR/Bd_contigs_virus_gene_hits.fasta -a $INDIR/Bd_contigs_virus.orfs.faa -m -o $INDIR/Bd_contigs_virus.orfs.txt
esl-sfetch --index $INDIR/Bd_contigs_virus.orfs.faa

#search
hmmsearch -E 1e-5 --domtbl $INDIR/Bd_virus.rep_gene.domtbl lib/profile/ssDNAvirus.rep.hmm $INDIR/Bd_contigs_virus.orfs.faa > $INDIR/Bd_virus.rep_gene.hmmer
ssearch36 -m 8c -E 1e-3 lib/viral_pep/Capsid.fa $INDIR/Bd_contigs_virus.orfs.faa > $INDIR/Bd_virus.Cap_gene.SSEARCH.tab

phmmer -E 1e-5 --domtbl $INDIR/Bd_virus.Cap_gene.domtbl lib/viral_pep/Capsid.fa $INDIR/Bd_contigs_virus.orfs.faa > $INDIR/Bd_virus.Cap_gene.phmmer

# extract
grep -v '^#' $INDIR/Bd_virus.rep_gene.domtbl | awk '{print $1}' | esl-sfetch -f $INDIR/Bd_contigs_virus.orfs.faa - > $INDIR/Bd_virus.rep_gene.fas

# $INDIR/rep_gene.refseq.fa
cat $INDIR/Bd_viral_hits.refseq_edit.rep_gene.fa  >> $INDIR/Bd_virus.rep_gene.fas

muscle -in $INDIR/Bd_virus.rep_gene.fas -out $INDIR/Bd_virus.rep_gene.afa

grep -v '^#' $INDIR/Bd_virus.Cap_gene.domtbl | awk '{print $1}' | esl-sfetch -f $INDIR/Bd_contigs_virus.orfs.faa - > $INDIR/Bd_virus.Cap_gene.fas

muscle -in $INDIR/Bd_virus.Cap_gene.fas -out $INDIR/Bd_virus.Cap_gene.afa

module load fasttree

FastTreeMP -gamma -lg < $INDIR/Bd_virus.rep_gene.afa > $INDIR/Bd_virus.rep_gene.ft.tre

FastTreeMP -gamma -lg < $INDIR/Bd_virus.Cap_gene.afa > $INDIR/Bd_virus.Cap_gene.ft.tre
