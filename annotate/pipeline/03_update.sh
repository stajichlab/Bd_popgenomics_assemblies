#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem 64gb
#SBATCH --time=3-00:00:00   
#SBATCH --output=logs/update.%a.log
#SBATCH --job-name="fun-update-Bd"
module unload perl
module unload miniconda2
module unload miniconda3
module load anaconda3
module load funannotate/development
module unload perl
module unload python
source activate funannotate

export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)

CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
 CPUS=2
fi

INDIR=genomes
OUTDIR=funannotate
mkdir -p $OUTDIR

SAMPFILE=samples.tab
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi
species_name=Batrachochytrium_dendrobatidis

sed -n ${N}p $SAMPFILE | while read SAMP
do
    name=$(echo "$SAMP" | perl -p -e 's/[\(\)]//g')
    PREFIX=$(echo "$name" | perl -p -e 's/\-/./g')
    if [ ! -f $INDIR/$name.masked.fasta ]; then
	echo "No genome for $INDIR/$name.masked.fasta yet - run 00_mash.sh $N"
	exit
    fi
    funannotate update -i $OUTDIR/$name --cpus $CPUS --sbt lib/Bd_pangenome.sbt
done
