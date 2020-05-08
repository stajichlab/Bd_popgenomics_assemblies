#!/usr/bin/bash
#SBATCH -p short -n 24 --mem 48gb --out logs/bwa_SE.%a.log -N 1

module load samtools/1.10
module load bwa
module load bam2fastq

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input_SE
OUT=aln
SAMPLEFILE=strains_SE.txt
TEMP=/scratch
#REF=genome/assembled_TF5a1.fa
REF=genome/Bd_genome_w_virus.fa

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $OUT

sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do
    ALNFILE=$OUT/$STRAIN.cram
    echo "processing $STRAIN"
    if [ ! -f $ALNFILE ]; then
	    pbzip2 -dc $INDIR/${LIBRARY}.fastq.bz2 > $TEMP/${LIBRARY}.fastq
	    bwa mem -t $CPU $REF /scratch/${LIBRARY}.fastq | samtools view -Obam - | samtools sort --threads $CPU --reference $REF -Ocram -o $ALNFILE -
	    rm $TEMP/${LIBRARY}.fastq
    fi
done
