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
 
  #print(ticks)
  pdffile= sprintf("plots/%s.chromplot.pdf",stem)
  Title = sprintf("%s PanGene Count",stem)
  d$Chromosome <- d$Chr
  #sc <- scale_colour_gradientn(colours = myPalette(10), limits=c(1, MAX_SCAFFOLDS))
  p <- ggplot(d, aes(START, GENOME_COUNT)) + geom_point(aes(color=factor(floor(GENOME_COUNT/100))),
                                                      alpha=0.5,size=2) +
    labs(title=Title,xlab="Chromosome") +
    facet_grid(Chromosome ~ .) +
#    scale_y_continuous(expand=c(0,0)) + 
    #    scale_colour_gradientn(colours = terrain.colors(5)) +
    #      scale_color_manual(values=rep(c("black","grey","black","grey"))) +
          scale_colour_brewer(palette = "Dark2") + 
       theme_cowplot(12)
  ggsave(pdffile,p,width=10,height=10)
}


