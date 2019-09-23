#!?usr/bin/bash
#SBATCH -p short --out make_tree.log

module load IQ-TREE

if [ ! -f virus_haplotypes.haploid.fasaln ]; then
	ln -s ../virus_haplotypes.haploid.fasaln .
fi
iqtree -s virus_haplotypes.haploid.fasaln -nt 1 -bb 1000 -alrt 1000
