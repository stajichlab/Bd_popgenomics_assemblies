#!/usr/bin/bash
#SBATCH --mem 64gb -p batch --out logs/metaspades.%a.log -n 8 -N 1

module load SPAdes
MEM=64

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=ref_aln
OUT=asm
mkdir -p $OUT
SAMPLEFILE=samples_PE.dat

sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN
do
	OUTNAME=$OUT/$(echo "$STRAIN" | perl -p -e 's/[\(\)]//g')
	if [ ! -d $OUTNAME ]; then
		metaspades.py -o $OUTNAME --mem $MEM -t $CPU -1 $INDIR/${STRAIN}_1.fastq.gz -2 $INDIR/${STRAIN}_2.fastq.gz
	elif [ ! -f $OUTNAME/scaffolds.fasta ]; then
		metaspades.py -o $OUTNAME --continue -t $CPU
	else
		echo "already run metaspades on $OUT/$STRAIN"
	fi
done
