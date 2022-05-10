#!/usr/bin/bash
#SBATCH -p batch -N 1 -n 32 --mem 128gb  --out pirate.log -J PIRATE

module unload miniconda3
module unload miniconda2
module unload python
module unload perl
module unload anaconda2
module load anaconda3
export TEMPDIR=/scratch
export TMPDIR=/scratch
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

source activate pirate

~/projects/PIRATE/bin/PIRATE -f mRNA -i gff3 -t $CPUS -n -a --rplots -o Bd_PIRATE_3
