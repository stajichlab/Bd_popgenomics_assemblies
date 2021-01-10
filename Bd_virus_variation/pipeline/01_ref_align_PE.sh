#!/usr/bin/bash
#SBATCH -p short -n 24 --mem 64gb --out logs/bwa_PE.%a.log -N 1

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

INDIR=input_PE
OUT=aln
SAMPLEFILE=strains_PE_all.txt
TEMP=/scratch
#REF=genome/assembled_TF5a1.fa
REF=genome/Bd_genome_w_virus.fa

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $OUT

sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do
	echo "$LIBRARY"
    ALNFILE=$OUT/$STRAIN.cram
    if [ ! -f $ALNFILE ]; then
	    bwa mem -t $CPU $REF $INDIR/${LIBRARY}_?.fastq.gz | samtools fixmate --threads $CPU -Ocram --reference $REF - $TEMP/${LIBRARY}.fix.cram
	    samtools sort --threads $CPU --reference $REF -Ocram -o $ALNFILE $TEMP/${LIBRARY}.fix.cram
	    rm $TEMP/${LIBRARY}.fix.cram
    fi
done
