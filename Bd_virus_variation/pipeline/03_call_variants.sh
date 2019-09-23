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
GENOME=genome/assembled_TF5a1.fa
for ploidy in haploid diploid
do
    RESULT=$OUTDIR/BdVirus.$ploidy.vcf.gz
    FILTERED=$OUTDIR/BdVirus.$ploidy.filtered.vcf.gz
    if [ $ploidy == "haploid" ]; then
	bcftools mpileup -Ou -f $GENOME $INDIR/*.bam | bcftools call --ploidy 1 -vmO z -o $RESULT
    else
	bcftools mpileup -Ou -f $GENOME $INDIR/*.bam | bcftools call -vmO z -o $RESULT
    fi
    tabix -p vcf $RESULT
    bcftools stats -F $GENOME -s - $RESULT > $RESULT.stats
    bcftools filter -O z -o $FILTERED -s LOWQUAL -s LOWQUAL -i'%QUAL>10' $RESULT
    tabix -p vcf $FILTERED
done
