pushd genomes
while read ID STRAIN; do 
	if [[ -f ../../asm_single/genomes/$ID.cleaned.fasta && ! -f $STRAIN.cleaned.fasta ]]; then  
		echo "ln -s ../../asm_single/genomes/$ID.cleaned.fasta $STRAIN.cleaned.fasta"; 
	fi; 
done < ../../asm_single/samples.all.dat
