#!/usr/bin/bash
#SBATCH -p short -n 24 --mem 16gb --out logs/bwa_refaln.%a.log -N 1

module load samtools/1.9
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

INDIR=reads
OUT=ref_aln
SAMPLEFILE=read_files_SE.dat
TEMP=/scratch
REF=ref_genome/Bd_JEL423.fasta
if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $OUT

sed -n ${N}p $SAMPLEFILE | while read STRAIN
do
    BAM=$OUT/$STRAIN.unmapped.bam
    if [ ! -f $BAM ]; then
	    pbzip2 -dc $INDIR/${STRAIN}.fastq.bz2 > $TEMP/${STRAIN}.fastq
	    bwa mem -t $CPU $REF /scratch/${STRAIN}.fastq | samtools view -h -u -f4 -Obam - | samtools sort --threads $CPU -Obam -o $BAM -
	    rm $TEMP/${STRAIN}.fastq
    fi
    if [ ! -f $OUT/${STRAIN}_1.fastq.gz ]; then
	bam2fastq --no-aligned -o $OUT/${STRAIN}#.fastq $BAM
	find $OUT -name "${STRAIN}_*" -size 0 | xargs rm
	pigz $OUT/${STRAIN}_[12M].fastq
    fi
done
