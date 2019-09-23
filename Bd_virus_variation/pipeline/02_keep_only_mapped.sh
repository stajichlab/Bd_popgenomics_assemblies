#!/usr/bin/bash
#SBATCH -N 1 -n 24 -p short --out logs/keep_only_mapped.log
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
module load samtools
module load parallel
module unload perl
mkdir -p mapped

#parallel -j $CPU samtools view --threads 2 -O bam -o mapped/{/.}.bam --reference genome/assembled_TF5a1.fa -F 4 {} ::: aln/*.cram
parallel -j $CPU samtools index {} ::: mapped/*.bam
