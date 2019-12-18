#!/usr/bin/bash

#SBATCH -p short -N 1 -n 2 --mem 8gb

perl scripts/summarize_annotation.pl > annotation_summary.tsv
