
covdat = read.csv("Bd_coverage_virus.tsv",header=T,sep="\t",row.names=1)
#summary(covdat)

pdf("Bd_MTCov.pdf")

colct <- ncol(covdat)

dropcols <- c("DS022322.1","TF5a1_Bdvirus")
MTchr <- covdat$DS022322.1
# remove MT and 
nuc <- covdat[,!(names(covdat) %in% dropcols)]
head(nuc)
nucavg <- rowSums(nuc)/ncol(nuc)
MTnucRatio <- MTchr / nucavg

hist(MTchr,50,main="MT coverage distribution")
hist(nucavg,50,main="Nuclear coverage distribution")
hist(MTnucRatio,50,main="MT:Nuclear ratio")

pdf("Bd_VirCov.pdf")
vironly <- subset(covdat,covdat$TF5a1_Bdvirus > 0)
nrow(vironly)

vir <- vironly$TF5a1_Bdvirus

nucVirOnly <- vironly[,!(names(vironly) %in% dropcols)]
nucVirAvg <- rowSums(nucVirOnly)/ncol(nucVirOnly)

DF <- data.frame(NucAvgCov = nucVirAvg,
                 VirCov    = vir,
                 VirRatio     = vir / nucVirAvg)

write.csv(DF,"virus_hits.csv")
rownames(DF) = rownames(vironly)

plot(VirCov ~ NucAvgCov,data=DF,main="Virus to Nuclear Coverage",pch=16,cex=0.1)
text(VirCov ~ NucAvgCov, labels=rownames(DF),data=DF,font=1,cex=0.7,
     main="Virus:Nuclear coverage")

plot(VirCov ~ NucAvgCov,data=DF,main="Virus to Nuclear Coverage",pch=2,cex=0.5)

boxplot(DF$VirRatio,outpch=NA,main="Virus:Nuclear Ratio")
stripchart(DF$VirRatio,
           vertical = TRUE, method = "jitter", 
           pch = 21, col = "maroon", bg = "bisque", 
           add = TRUE) 





