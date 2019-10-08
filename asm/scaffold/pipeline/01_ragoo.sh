#!/usr/bin/bash
#SBATCH -p short  --mem 16gb -N 1 -n 8 --out logs/ragooo.%a.log

module unload miniconda3
module load anaconda3
module load minimap2

module load AAFTF
source activate ragoo

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -z $CPU ]; then
	CPU=1
fi
SAMPLEFILE=samples.dat
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
NAME=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
ASM=genomes
OUT=scaffolds
CTGS=$(realpath $ASM/$NAME.sorted_shovill.fasta)
REF=$(realpath Public_assemblies/genome/GCA_000149865.1_BD_JEL423_genomic.fna)
echo $REF
echo $CTGS
SCAF=$OUT/$NAME

LEFT=$WORKDIR/${BASE}_filtered_1.fastq.gz
RIGHT=$WORKDIR/${BASE}_filtered_2.fastq.gz

mkdir -p $SCAF
pushd $SCAF
ln -s $CTGS contigs.fasta
ln -s $REF reference.fasta
ragoo.py -m $(which minimap2) -t $CPU -b -C contigs.fasta reference.fasta
AAFTF assess -i ragoo_output/ragoo.fasta -r assembly.stats.txt
popd
