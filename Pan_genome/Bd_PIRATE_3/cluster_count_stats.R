library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(cowplot)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# read in the cluster table
topdir = "per-strain-pangenome"
file.names = list.files(path = topdir, pattern = "*.pangene_depth.tab")
MAX_SCAFFOLDS=20
for (filei in 1:length(file.names)) {
  file = file.names[filei]
  stem = sub("Batrachochytrium_dendrobatidis_","",file)
  stem = sub(".pangene_depth.tab","",stem)
  df <- read.csv(sprintf("%s/%s",topdir,file),header=T,sep="\t",stringsAsFactors=FALSE)
  df$Chr <- sub("scaffold_","",df$CHROM,perl=TRUE)
#  chrlist = c(1:8)
  chrlist = c(1:MAX_SCAFFOLDS)
  d = df[df$Chr %in% chrlist, ]
  d <- d[order(as.numeric(d$Chr), as.numeric(d$START)), ]
  d$GENOME_COUNT = as.numeric(d$GENOME_COUNT)
  d$START = as.numeric(d$START)
  d$Chr = as.numeric(d$Chr)
  d$index = rep.int(seq_along(unique(d$Chr)), 
                    times = tapply(d$START,d$Chr,length))

  d$pos=NA
  nchr = length(unique(d$Chr))
  lastbase=0
  ticks = NULL
  minor = vector(,8)
  for (i in chrlist ) {
   if (i==1) {
      d[d$index==i, ]$pos = d[d$index==i, ]$START
    } else {
      ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
      lastbase = lastbase + max(d[d$index==(i-1),"START"])
      minor[i] = lastbase
      d[d$index == i,"START"] =
          d[d$index == i,"START"] - min(d[d$index==i,"START"]) + 1
      d[d$index == i,"END"] = lastbase
      d[d$index == i, "pos"] = d[d$index == i,"START"] + lastbase
    }
  }

  ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
  print(ticks)
  pdffile= sprintf("plots/Batrachochytrium_dendrobatidis_%s.chromplot.pdf",stem)
  Title = sprintf("%s PanGene Count",stem)
  d$Chromosome <- d$Chr
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, MAX_SCAFFOLDS))
  p <- ggplot(d, aes(pos, GENOME_COUNT)) + geom_point(aes(color=Chromosome),
                                               alpha=0.5,size=0.5) +
    labs(title=Title,xlab="Chromosome",scales="free_y") +
    scale_x_continuous(name="Chromosome", expand = c(0, 0),
                       breaks=ticks,
                       labels=(unique(d$Chr))) +
      scale_y_continuous(expand=c(0,0)) + 
#    scale_colour_gradientn(colours = terrain.colors(5)) +
#      scale_color_manual(values=rep(c("black","grey","black","grey"))) +
#      scale_colour_brewer(palette = "Dark2") + 
    theme(legend.position="none", 
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + 
    theme_cowplot(12)
  ggsave(pdffile,p,width=7,height=3)
  break
}


