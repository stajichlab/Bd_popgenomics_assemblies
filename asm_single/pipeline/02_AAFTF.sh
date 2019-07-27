#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 48gb -J BdAAFTF --out logs/AAFTF.%a.log -p intel --time 64:00:00
source ~/.bashrc
hostname
MEM=48
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
module load AAFTF

INPUT=input
SAMPLEFILE=samples.dat
PREFIX=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
ASM=genomes
TMPDIR=/scratch/$USER
MINLEN=500

mkdir -p $ASM

if [ -z $CPU ]; then
    CPU=1
fi

ASMFILE=$ASM/${BASE}.spades.fasta
WORKDIR=working_AAFTF
VECCLEAN=$ASM/${BASE}.vecscreen.fasta
PURGE=$ASM/${BASE}.sourpurge.fasta
CLEANDUP=$ASM/${BASE}.rmdup.fasta
PILON=$ASM/${BASE}.pilon.fasta
SORTED=$ASM/${BASE}.sorted.fasta
STATS=$ASM/${BASE}.sorted.stats.txt
STATS=$ASM/${BASE}.vecclean.stats.txt

#LEFT=$INPUT/${PREFIX}_1.fastq.gz
#RIGHT=$INPUT/${PREFIX}_2.fastq.gz
LEFT=$WORKDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$WORKDIR/${PREFIX}_filtered_2.fastq.gz

mkdir -p $WORKDIR

echo "$BASE"
if [ ! -f $ASMFILE ]; then    
    if [ ! -f $LEFT ]; then
	echo "Cannot find LEFT $LEFT or RIGHT $RIGHT - did you run"
	echo "$OUTDIR/${BASE}_R1.fq.gz $OUTDIR/${BASE}_R2.fq.gz"
	exit
    fi
    spades.py -s $LEFT --threads $CPU --mem $MEM -o  $WORKDIR/spades_${BASE} 

    if [ -f $WORKDIR/spades_${BASE}/scaffolds.fa ]; then
	rsync -av $WORKDIR/spades_${BASE}/scaffolds.fa $ASMFILE
    else	
	echo "Cannot find $OUTDIR/spades_${BASE}/scaffolds.fa"
    fi
    
    if [ -s $ASMFILE ]; then
	rm -rf $WORKDIR/spades_${BASE}
    else
	echo "SPADES must have failed, exiting"
	exit
    fi
fi

if [ ! -f $VECCLEAN ]; then
    AAFTF vecscreen -i $ASMFILE -c $CPU -o $VECCLEAN 
fi
if [ ! -f $STATS ]; then
	AAFTF assess -i $VECCLEAN -r $STATS
fi
exit
