cd genomes
while read ID STRAIN; do 
	FIXNAME=$(echo "$STRAIN" | perl -p -e 's/[\(\)]//g')

	if [[ -f ../../asm/genomes/$ID.cleaned.fasta && ! -f $FIXNAMAE.cleaned.fasta ]]; then  
		echo "ln -s ../../asm/genomes/$ID.cleaned.fasta $STRAIN.cleaned.fasta"; 
	fi; 
done < ../../asm/samples.dat
cd ..


