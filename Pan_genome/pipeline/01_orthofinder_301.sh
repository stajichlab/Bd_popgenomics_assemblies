#!/usr/bin/bash
#SBATCH -p intel --mem 64gb -N 1 -n 32 --out logs/orthofinder.301.%A.diamond.log -J orthofinder --time 72:00:00

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
orthofinder.py -f $DIR -a $CPUS -t $CPUS -S diamond -o $OUTDIR/OrthoFinder.Bd_301 -n Bd301 -og
