#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 16 --nodes 1 --mem 24G --out logs/predict.%a.log
conda activate
module unload miniconda2
module unload miniconda3
module load funannotate/development
module unload perl
module unload python
source activate funannotate
#export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)

which python
GMFOLDER=`dirname $(which gmhmme3)`
#genemark key is needed
if [ ! -f ~/.gm_key ]; then
	ln -s $GMFOLDER/.gm_key ~/.gm_key
fi

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
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
species="Batrachochytrium dendrobatidis"
SEED_SPECIES="batrachochytrium_dendrobatidis_jel423"
# use larger datasets for prediction?
TRANSCRIPTS=$(realpath lib/Trinity_combined.fasta)
#Trinity_GG_SpoZoo_PE.fasta
augustus_model=batrachochytrium_dendrobatidis_jel423
PEPS=$(realpath lib/informant.aa)
sed -n ${N}p $SAMPFILE | while read SAMP
do
	name=$(echo "$SAMP" | perl -p -e 's/[\(\)]//g')
	PREFIX=$(echo "$name" | perl -p -e 's/\-/./g')
	if [ ! -f $INDIR/$name.masked.fasta ]; then
		echo "No genome for $INDIR/$name.masked.fasta yet - run 00_mash.sh $N"
		exit
	fi
	echo "$name is name"
 	mkdir $name.predict.$$
 	pushd $name.predict.$$
    	funannotate predict --cpus $CPU --keep_no_stops --SeqCenter NCBI --busco_db fungi_odb9 --strain "$SAMP" \
      -i ../$INDIR/$name.masked.fasta --protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
      -s "$species"  -o ../$OUTDIR/$name --busco_seed_species $SEED_SPECIES --name="$PREFIX" \
	 --min_training_models 50 --transcript_evidence $TRANSCRIPTS \
      --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH --augustus_species $augustus_model

	popd
 	rmdir $name.predict.$$
done
