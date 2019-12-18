#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 64G -p intel --out logs/train.%a.log -J trainFun --time 48:00:00
conda activate
module unload miniconda2
module unload miniconda3
module load funannotate/development
module unload python
module unload perl
source activate funannotate
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
export PASAHOME=`dirname $(which Launch_PASA_pipeline.pl)`

CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPUS=2
fi

# use the smaller paired-end dataset for training
RNASEQ=lib/Trinity_GG_SpoZoo_PE.fasta
LEFT=rnaseq/Bd_paired_R1.fq.gz
RIGHT=rnaseq/Bd_paired_R2.fq.gz

INDIR=genomes
OUTDIR=funannotate
SAMPFILE=samples.tab
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $(expr $MAX) ]; then
    MAXSMALL=$(expr $MAX)
    echo "$N is too big, only $MAXSMALL lines in $SAMPFILE"
    exit
fi
INTRONLEN=3000
SPECIES="Batrachochytrium dendrobatidis"

sed -n ${N}p $SAMPFILE | while read SAMP
do
	name=$(echo "$SAMP" | perl -p -e 's/[\(\)]//g')
	ISOLATE=$(echo $name | perl -p -e 's/\-/./g')

	SORTED=$INDIR/$name.masked.fasta
	if [ ! -s $SORTED ]; then
		echo "Cannot find $SORTED"
		exit
	fi
	funannotate train -i $SORTED --species "$SPECIES" --isolate $ISOLATE --cpus $CPUS \
    	-o $OUTDIR/$ISOLATE --max_intronlen $INTRONLEN --trinity $RNASEQ --jaccard_clip --left $LEFT --right $RIGHT 
done
