library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(cowplot)
# read in the cluster table

df <- read.csv("per-strain-pangenome/Batrachochytrium_dendrobatidis_UM142.pangene_depth.tab",
                 header=T,sep="\t")
df$Chr <- sub("scaffold_","",df$CHROM,perl=TRUE)
chrlist = c(1:8)
d=df[df$Chr %in% chrlist, ]
d <- d[order(d$Chr, d$START), ]
d$index = rep.int(seq_along(unique(d$Chr)), times = tapply(d$START,d$Chr,length))

d$pos=NA

nchr = length(unique(d$Chr))
lastbase=0
ticks = NULL
minor = vector(,8)
for (i in 1:8 ) {
  if (i==1) {
    d[d$index==i, ]$pos=d[d$index==i, ]$START
  } else {
    ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
    lastbase = lastbase + max(d[d$index==(i-1),"START"])
    minor[i] = lastbase
    d[d$index == i,"START"] =
        d[d$index == i,"START"]-min(d[d$index==i,"START"]) +1
    d[d$index == i,"END"] = lastbase
    d[d$index == i, "pos"] = d[d$index == i,"START"] + lastbase
  }
}

ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
print(ticks)
#pdffile= 'plots/Batrachochytrium_dendrobatidis_UM142.pdf'
Title = "PanGene Count"
d$Chromosome <- d$Chr

p <- ggplot(d, aes(pos, GENOME_COUNT)) + geom_point(aes(color=Chromosome),
                                               alpha=0.5,size=0.5) +
  labs(title=Title,xlab="Chromosome",scales="free_y") +
  scale_x_continuous(name="Chromosome", expand = c(0, 0),
                     breaks=ticks,
                     labels=(unique(d$Chr))) +
    scale_y_continuous(expand=c(0,0))+
    scale_colour_brewer(palette = "Dark2") + theme_minimal() +
   theme(legend.position="none", panel.border = element_blank(),panel.grid.minor = element_blank(),
         panel.grid.major = element_blank()) + theme_cowplot(12)
p
#ggsave(pdffile,p,width=7,height=3)


