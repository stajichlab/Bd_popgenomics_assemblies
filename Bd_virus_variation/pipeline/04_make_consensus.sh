#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --out logs/consensus.log

module load samtools
module load bcftools

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

VCFDIR=vcf
HAPLOTYPES=virus_haplotypes
mkdir -p $HAPLOTYPES
GENOME=genome/Bd_genome_w_virus.fa

for ploidy in haploid diploid
do
    FILTERED=$VCFDIR/BdVirus.$ploidy.filtered.vcf.gz
    bcftools norm -f $GENOME -c s --threads $CPUS -Oz -o $VCFDIR/BdVirus.$ploidy.norm.vcf.gz $FILTERED
    bcftools index -f $VCFDIR/BdVirus.$ploidy.norm.vcf.gz
    for strain in $(bcftools query -l $VCFDIR/BdVirus.$ploidy.norm.vcf.gz);
    do
		   strainout=$(basename $strain .bam)
       echo "$strainout"
	     bcftools consensus -f $GENOME --haplotype A --sample $strain $VCFDIR/BdVirus.$ploidy.norm.vcf.gz | perl -p -e "s/>/>$strainout /" > $HAPLOTYPES/$strainout.$ploidy.fa
     done
     cat virus_haplotypes/*.$ploidy.fa > virus_haplotypes.$ploidy.fa
     muscle -clw -in virus_haplotypes.$ploidy.fa -out virus_haplotypes.$ploidy.aln -quiet
     muscle -in virus_haplotypes.$ploidy.fa -out virus_haplotypes.$ploidy.fasaln -quiet

done
