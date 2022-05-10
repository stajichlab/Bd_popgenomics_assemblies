#!/usr/bin/bash
#SBATCH -p short --mem 64gb -N 1 -n 32 --out logs/orthofinder.mmseq.log -J orthofinder

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
DIR=test
JOBS=$(expr $CPUS / 2)
OUTDIR=results
OUTEXT=test
N=1
orthofinder.py -f test -t $CPUS -S mmseqs -o $OUTDIR/OrthoFinder.Run$N -n $OUTEXT
