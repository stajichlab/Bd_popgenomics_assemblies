#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 48gb -J shovill --out logs/shovill.%a.log -p intel --time 64:00:00
source ~/.bashrc
module load AAFTF
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


INPUT=input
SAMPLEFILE=samples.dat
PHYLUM=Chytridiomycota
PREFIX=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
ASM=genomes
TMPDIR=/scratch/$USER
MINLEN=500

mkdir -p $ASM

if [ -z $CPU ]; then
    CPU=1
fi

ASMFILE=$ASM/${BASE}.spades_shovill.fasta
WORKDIR=working_AAFTF
VECCLEAN=$ASM/${BASE}.vecscreen_shovill.fasta
PURGE=$ASM/${BASE}.sourpurge_shovill.fasta
CLEANDUP=$ASM/${BASE}.rmdup_shovill.fasta
PILON=$ASM/${BASE}.pilon_shovill.fasta
SORTED=$ASM/${BASE}.sorted_shovill.fasta
STATS2=$ASM/${BASE}.sorted_shovill.stats.txt
STATS=$ASM/${BASE}.vecclean_shovill.stats.txt

LEFTORIG=$INPUT/${PREFIX}_1.fastq.gz
RIGHTORIG=$INPUT/${PREFIX}_2.fastq.gz
LEFT=$WORKDIR/${PREFIX}_filtered_1.fastq.gz
RIGHT=$WORKDIR/${PREFIX}_filtered_2.fastq.gz

mkdir -p $WORKDIR

echo "$BASE"
if [[ ! -f $ASMFILE && ! -f $ASMFILE.bz2 ]]; then    
    if [ ! -f $LEFT ]; then
	echo "Cannot find LEFT $LEFT or RIGHT $RIGHT - did you run"
	echo "$OUTDIR/${BASE}_R1.fq.gz $OUTDIR/${BASE}_R2.fq.gz"
	exit
    fi
    module unload miniconda2
    module load miniconda3
    module unload perl
    source activate shovill

    shovill --cpu $CPU --ram $MEM --outdir $WORKDIR/shovill_${BASE} \
	--R1 $LEFT --R2 $RIGHT --depth 90 --tmpdir $TMPDIR --minlen $MINLEN --nocorr

    if [ -f $WORKDIR/shovill_${BASE}/contigs.fa ]; then
	rsync -av $WORKDIR/shovill_${BASE}/contigs.fa $ASMFILE
    else	
	echo "Cannot find $WORKDIR/shovill_${BASE}/contigs.fa"
    fi
    conda deactivate 
    
    if [ -s $ASMFILE ]; then
	rm -rf $WORKDIR/shovill_${BASE}
    else
	echo "SPADES must have failed, exiting"
	exit
    fi
fi

if [ ! -f $VECCLEAN ]; then
	if [[ -f $ASMFILE.bz2 && ! -f $ASMFILE ]]; then
		pbzip2 -k $ASMFILE.bz2
	fi
    AAFTF vecscreen -i $ASMFILE -c $CPU -o $VECCLEAN 
fi
if [ ! -f $STATS ]; then
	AAFTF assess -i $VECCLEAN -r $STATS
fi

if [ ! -f $PURGE ]; then
    AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $LEFTORIG  --right $RIGHTORIG
fi

if [ ! -f $CLEANDUP ]; then
   AAFTF rmdup -i $PURGE -o $CLEANDUP -c $CPU -m 500
fi

if [ ! -f $PILON ]; then
   AAFTF pilon -i $CLEANDUP -o $PILON -c $CPU --left $LEFTORIG  --right $RIGHTORIG
fi

if [ ! -f $PILON ]; then
    echo "Error running Pilon, did not create file. Exiting"
    exit
fi

if [ ! -f $SORTED ]; then
    AAFTF sort -i $PILON -o $SORTED
fi

if [ ! -f $STATS2 ]; then
    AAFTF assess -i $SORTED -r $STATS2
fi
