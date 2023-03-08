#!/usr/bin/env Rscript

# Plot Toxin Loci
# Author: Pedro G. Nachtigall

### Load and install Packages ###
packages<-c("ggplot2","gggenes", "ggrepel", "tidyverse","argparse")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
 install.packages(packages[!installed_packages])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

### SETUP ARGPARSE ###
parser <- ArgumentParser()
parser$add_argument("-i","--input", default="toxin_annotation_GENE.tsv", 
                    help="toxin annotation file with gene annotation as output by fromCDStoGENE.py [default: \"%(default)s\"]")
parser$add_argument("-o","--output", default="output",
                    help="output prefix to use [default: \"%(default)s\"], the default will output the file \"output_toxin_loci.pdf\" ")
args <- parser$parse_args()

### Plot Function###
ToxinPlot <- function(df=df,print=TRUE){
  A<-ggplot(df, aes(xmin = start, xmax = end, y = molecule, fill = family, forward = direction)) +
  		geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  		ggrepel::geom_text_repel(data = df, aes(x = start, y = molecule, label = gene), inherit.aes = F, nudge_y = 1) +
  		facet_wrap(~ molecule, scales = "free", ncol = 1) +
  		ylab('') +
  		xlab('') +
  		theme_genes()
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}


### READ DATA ###
df <- read.table(normalizePath(args$input), header=FALSE, sep = "\t")
colnames(df) <- c("gene","family","molecule","start","end","strand")
df <- df %>% mutate(strand = ifelse(strand == '+','forward','reverse')) %>% mutate(direction = ifelse(strand == 'forward',1,-1))

### MAKE PLOTS ###
pdf(paste0(args$output,"_toxin_loci.pdf"), width = 12, height=10)
ToxinPlot(df=df)
dev.off()

# END
