#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 8 --nodes 1 --mem 12G 
#SBATCH --out logs/mask.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=genomes
LIBRARY=lib/Bd_rmlib.fasta
SAMPFILE=samples.dat
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


sed -n ${N}p $SAMPFILE | while read SRA name
do
    GENOMEFILE=$INDIR/${name}.cleaned.fasta
    if [ ! -f $GENOMEFILE ]; then
	echo "Cannot find $name (${name}.cleaned.fasta) in $INDIR - may not have been run yet"
	exit
    fi
    
    if [ ! -f $OUTDIR/${name}.masked.fasta ]; then	
	module load funannotate/git-live
	module unload rmblastn
	module load ncbi-rmblast/2.9.0-p2
	export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
        
	if [ -f repeat_library/${name}.repeatmodeler-library.fasta ]; then
	    LIBRARY=repeat_library/${name}.repeatmodeler-library.fasta
	fi
	LIBRARY=$(realpath $LIBRARY)
	mkdir $name.mask.$$
	pushd $name.mask.$$
	funannotate mask --cpus $CPU -i ../$GENOMEFILE -o ../$OUTDIR/${name}.masked.fasta -l $LIBRARY
	mv funannotate-mask.log ../logs/${name}.mask.log
	popd
	rmdir $name.mask.$$
    else 
     echo "Skipping ${name} as masked already"
    fi
done
