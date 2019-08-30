#!/usr/bin/bash
#SBATCH -N 1 -n 16 -p short --mem 16gb --out logs/make_cov.%a.log

module load bwa
module load samtools/1.9
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi


if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then
 echo "need to provide a number by --array or cmdline"
 exit
fi
mkdir -p bam
SAMPLEFILE=samples.dat
PREFIX=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f2)

ASSEMBLY=genomes/$BASE.vecscreen.fasta
if [ ! -f $ASSEMBLY ]; then
    echo "No $ASSEMBLY file"
    exit
fi
BAM=bam/$BASE.remap.bam
FWD=input/${PREFIX}.fastq.gz

if [ ! -f $BAM ]; then
    if [ ! -f $ASSEMBLY.bwt ]; then
	bwa index $ASSEMBLY
    fi
    bwa mem -t $CPU $ASSEMBLY $FWD | samtools sort --threads $CPU -T /scratch -O bam -o $BAM -
fi

if [ ! -f $BAM.bai ]; then
    samtools index $BAM
fi
