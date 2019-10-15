#!/usr/bin/bash

VDIR=virus_search
mkdir -p $VDIR

for f in $(cat lib/ncbi_refseq_virus.txt)
do
	b=$(basename $f .gz)
	if [ ! -f $VDIR/$b ]; then
		curl -o $VDIR/$f ftp://ftp.ncbi.nih.gov/refseq/release/viral/$f
		pigz -d $VDIR/$f
	fi
done

