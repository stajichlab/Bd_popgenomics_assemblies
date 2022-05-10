#!/usr/bin/bash
#SBATCH -p short --mem 64gb -N 1 -n 64 --out logs/orthofinder.299.%A.diamond.log -J orthofinder

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

module unload miniconda3
module load miniconda2
module load orthofinder

# if we are on short/intel queue 
# need to figure out how to better determine this
module unload MMseqs2/10-6d92c
module load MMseqs2/10-6d92c-avx2
DIR=input
JOBS=$(expr $CPUS / 2)
OUTDIR=results
N=diamond
orthofinder.py -f $DIR -a $CPUS -t $CPUS -S diamond -o $OUTDIR/OrthoFinder.Bd_299 -n Bd299 -og
