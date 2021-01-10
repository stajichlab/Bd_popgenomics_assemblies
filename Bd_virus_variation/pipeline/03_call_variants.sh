#!/usr/bin/bash
#SBATCH -N 1 -n 8 -p short --out logs/call_variants.log
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
module load samtools
module load bcftools

mkdir -p vcf
OUTDIR=vcf
INDIR=mapped
GENOME=genome/Bd_genome_w_virus.fa
VIRHITS=virus_hits.csv
BAMFILES=virus_bamfiles.txt
SAMPLEFILES=virus_samples.txt
PLOIDYFILE=virus_haploid.txt
cut -d, -f1 $VIRHITS | tail -n +2 | perl -p -e 's/"//g; s/(\S+)/mapped\/$1.bam/' > $BAMFILES
cut -d, -f1 $VIRHITS | tail -n +2 | perl -p -e 's/"//g; s/(\S+)/mapped\/$1.bam\t$1/' > $SAMPLEFILES
cut -d, -f1 $VIRHITS | tail -n +2 | perl -p -e 's/"//g; s/(\S+)/$1\t1/' > $PLOIDYFILE
for ploidy in haploid diploid
do
    RESULT=$OUTDIR/BdVirus.$ploidy.vcf.gz
    FILTERED=$OUTDIR/BdVirus.$ploidy.filtered.vcf.gz
    if [ $ploidy == "haploid" ]; then
			bcftools mpileup --threads $CPU -Ou -f $GENOME --bam-list virus_bamfiles.txt -S $SAMPLEFILES | bcftools call -S $PLOIDYFILE -vmO z -o $RESULT
    else
			bcftools mpileup --threads $CPU -Ou -f $GENOME --bam-list virus_bamfiles.txt  -S $SAMPLEFILES | bcftools call -vmO z -o $RESULT
    fi
    tabix -p vcf $RESULT
    bcftools stats -F $GENOME -s - $RESULT > $RESULT.stats
    bcftools filter -O z -o $FILTERED -s LOWQUAL -s LOWQUAL -i'%QUAL>10' $RESULT
    tabix -p vcf $FILTERED
done
