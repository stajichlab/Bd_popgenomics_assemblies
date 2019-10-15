#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 2gb --out logs/remove_viral_contigs.%a.log

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi


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
QUERY=lib/assembled_TF5a1.fa
SORTED=$ASM/${BASE}.sorted_shovill.fasta
CLEAN=$ASM/${BASE}.cleaned.fasta
REPORTDIR=viral_contig_search
REPORTFILE=$REPORTDIR/${BASE}.viralSearch.FASTA.m9
mkdir -p $REPORTDIR

if [ ! -f $SORTED ]; then
    echo "No $SORTED file, need to run 02_AAFTF_shovill.sh $N"
	exit
fi

module load hmmer
module load fasta

if [ ! -f $REPORTFILE  ]; then
    fasta36 -E 1e-20 -T $CPU -m 8 $QUERY $SORTED > $REPORTFILE
fi

if [ -s $REPORTFILE ]; then
	if [[ ! -f $SORTED.ssi || $SORTED -nt $SORTED.ssi ]]; then
		esl-sfetch --index $SORTED
	fi
	CTGS=$REPORTDIR/${BASE}.viralSearch.contigs.fasta
	if [[ ! -s $CTGS || $REPORTFILE -nt $CTGS ]]; then
		cut -f2 $REPORTFILE | sort | uniq | esl-sfetch -f $SORTED - | perl -p -e "s/>/>${BASE}_/" > $CTGS
	fi
fi

if [[ ! -f $CLEAN || $REPORTFILE -nt $CLEAN ]]; then
	./scripts/remove_viral_contigs.py $REPORTFILE $SORTED > $CLEAN
fi
