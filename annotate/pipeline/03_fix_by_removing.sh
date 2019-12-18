#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem 8gb -p short
#SBATCH --time=2:00:00   
#SBATCH --output=logs/update_fix.%a.log
#SBATCH --job-name="updatefix-Bd"
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
	mkdir $name.fix.$$
 	pushd $name.fix.$$
	grep -v '^#' ../$OUTDIR/$name/update_results/*.models-need-fixing.txt | awk '{print $1}' > remove.${SLURM_ARRAY_TASK_ID}
	if [ -s remove.${SLURM_ARRAY_TASK_ID} ]; then
		funannotate fix -i ../$OUTDIR/$name/update_results/${species_name}_${name}.gbk -t  ../$OUTDIR/$name/update_results/${species_name}_${name}.tbl -d remove.${SLURM_ARRAY_TASK_ID}
	fi
	popd
	rm -rf $name.fix.$$
done
