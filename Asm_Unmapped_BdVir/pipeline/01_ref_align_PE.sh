#!/usr/bin/bash
#SBATCH -p short -n 24 --mem 16gb --out logs/bwa_refaln.%a.log -N 1

module load samtools/1.9
module load bwa
module load bam2fastq

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}
TEMP=/scratch
if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
OUT=ref_aln
SAMPLEFILE=samples_PE.dat
REF=ref_genome/Bd_JEL423.fasta
if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN
do
    BAM=$OUT/$STRAIN.unmapped.bam
    if [ ! -f $BAM ]; then
	bwa mem -t $CPU $REF $INDIR/${BASE}_[12].fastq.gz | samtools view -h -f8 -f4 -Obam - | samtools sort --threads $CPU -Obam -o $BAM -
    fi
    if [ ! -f $OUT/${STRAIN}_1.fastq.gz ]; then
	bam2fastq --no-aligned -o $OUT/${STRAIN}#.fastq $BAM
	find $OUT -name "${STRAIN}_*" -size 0 | xargs rm
	pigz $OUT/${STRAIN}_[12].fastq
    fi
done
