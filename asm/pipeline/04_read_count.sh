#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 24 --mem 24G -p short -J readcov --out logs/bbcount.%a.log --time 2:00:00
module load BBMap
hostname
MEM=24
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
SAMPLEFILE=samples.dat
ID=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
ASM=genomes
OUTDIR=mapping_report
SORTED=$(realpath $ASM/${BASE}.cleaned.fasta)

LEFT=$(realpath $INDIR/${ID}_1.fastq.gz)
RIGHT=$(realpath $INDIR/${ID}_2.fastq.gz)
mkdir -p $OUTDIR
BASE=$(echo $BASE | perl -p -e 's/[\)\(]//g')
if [ ! -s $OUTDIR/${BASE}.bbmap_covstats.txt ]; then
	mkdir -p N$N.$$.bbmap
	pushd N$N.$$.bbmap
	bbmap.sh -Xmx${MEM}g ref=\"$SORTED\" in=$LEFT in2=$RIGHT covstats=../$OUTDIR/${BASE}.bbmap_covstats.txt  statsfile=../$OUTDIR/${BASE}.bbmap_summary.txt
	popd
	rm -rf N$N.$$.bbmap
fi
