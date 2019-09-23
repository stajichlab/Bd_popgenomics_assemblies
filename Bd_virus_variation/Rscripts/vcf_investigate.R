library(vcfR)

vcf <- read.vcfR("vcf/BdVirus.diploid.filtered.vcf.gz")
dna_file <- "genome/assembled_TF5a1.fa"
gff_file <- "genome/assembled_TF5a1.gff"
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")
chrom <- create.chromR(name='TF5a1_Bdvirus', vcf=vcf, seq=dna, ann=gff)
pdf("TF5a1_Bdvirus.vcfR.pdf")
plot(chrom)
chromoqc(chrom, dp.alpha=20)
