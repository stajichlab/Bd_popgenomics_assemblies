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
SAMPLEFILE=reads_SE.dat

sed -n ${N}p $SAMPLEFILE | while read STRAIN
do
	if [ ! -d $OUT/$STRAIN ]; then
		spades.py -o $OUT/$STRAIN --mem $MEM -t $CPU -s $INDIR/${STRAIN}_M.fastq.gz
	elif [ ! -f $OUT/$STRAIN/scaffolds.fasta ]; then
		spades.py -o $OUT/$STRAIN --continue -t $CPU
	else
		echo "already run spades on $OUT/$STRAIN"
	fi
done
