#!/usr/bin/bash
#SBATCH -N 1 -n 24 -p short --out logs/keep_only_mapped.log
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
module load samtools/1.10
module load parallel
module unload perl
#GENOME=genome/assembled_TF5a1.fa
GENOME=genome/Bd_genome_w_virus.fa
CHROM=genome/Bd_genome_w_virus.chroms.bed
COV=coverage
MAP=mapped
IN=aln
mkdir -p $MAP $COV
parallel -j $CPU [[ ! -f {}.crai ]] \&\& samtools index {} ::: $IN/*.cram
parallel -j $CPU [[ ! -f $MAP/{/.}.bam ]] \&\& samtools view --threads 2 -O bam -o $MAP/{/.}.bam --reference $GENOME -F 4 {} TF5a1_Bdvirus ::: $IN/*.cram
parallel -j $CPU [[ ! -f {}.bai ]] \&\& samtools index {} ::: $MAP/*.bam
parallel -j $CPU [[ ! -f $COV/{.}.regions.bed.gz ]] \&\& mosdepth -f $GENOME -x -n --by $CHROM $COV/{/.} {} ::: $IN/*.cram

