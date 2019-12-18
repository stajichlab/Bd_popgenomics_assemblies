#!/bin/bash
#SBATCH -p short logs/find_missing_masked.log

CPU=1

INDIR=genomes
OUTDIR=genomes
SAMPFILE=samples.tab
IFS=,
N=1
mkdir -p empty
Species="Batrachochytrium_dendrobatidis"
echo "#!/usr/bin/bash" > move_$$.sh

m=$(cat $SAMPFILE | while read Strain
    do
	species=$(echo "$Species" | perl -p -e "chomp; s/$Strain//; s/\s+/_/g;")
	strain=$(echo $Strain | perl -p -e 'chomp; s/\s+/_/g')
	if [ ! -z $strain ]; then
 	    outname="${species}_$strain"
	else
	    outname="${species}"
	fi
	name=$Strain
	proteins=funannotate/${name}/annotate_results/$outname.proteins.fa
	updatepep=funannotate/${name}/update_results/$outname.proteins.fa
	if [ ! -f $INDIR/${name}.scaffolds.fasta ]; then
	    echo -e "\tCannot find ${name}.scaffolds.fasta in $INDIR - may not have been run yet ($N)" 1>&2
	elif [ ! -f $OUTDIR/${name}.masked.fasta ]; then
	    echo "need to run mask on $name ($N)" 1>&2
	elif [ ! -f $proteins ]; then
            echo "need to run annotate on $name ($N) no $proteins" 1>&2
	    echo $N
	elif [[ $updatepep -nt $proteins ]]; then
	    echo "need to run annotate on $name ($N) $updatepep is newer than $proteins" 1>&2
	    echo "mkdir -p funannotate/${name}/archive_$$" >> move_$$.sh
	    echo "mv funannotate/${name}/annotate_misc funannotate/${name}/annotate_results funannotate/${name}/archive_$$" >> move_$$.sh
	    echo $N
	fi
	
	N=$(expr $N + 1)
    done | uniq | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')

echo "sbatch --array=$m pipeline/04_function.sh"
