library(tidyverse)
library(DESeq2)
library(cowplot)

getwd()
sampleInfo = read.delim("shortRNAseq.txt")
sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)

## I am assuming feature counts finished
countdata = read.table("fly_counts.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
# or do it the tidy way
countdata = countdata[ ,6:ncol(countdata)]
# Remove crap from colnames in countdata

# format counts matrix colnames
temp = colnames(countdata)
temp = gsub("X.dfs6.pub.itamburi.ee283.week3.prob5.RNAbam.","",temp)
temp = gsub(".sort.bam","",temp)
colnames(countdata) = temp
# check everything matches
cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)


# create DEseq2 object & run DEseq
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds, contrast = c("TissueCode","B","E") )
summary(res)

# check dispersions and variance
plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )

###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)
sampleDists = dist( t( assay(rld) ) )

# heat map
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

# PCs
# wow you can sure tell tissue apart
print( plotPCA( rld, intgroup = "TissueCode") )

# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))



lfc = data.frame(res) %>%
  rownames_to_column(var = "gene.id") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(sig.label = ifelse(padj <=0.01, gene.id, NA))


ggplot(lfc, aes(log2FoldChange, -log10(padj))) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", alpha = .5) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", alpha = .5) +
  geom_point(size = .2, alpha = .4, color = "grey80") +
  geom_point(size = .2,
             data = subset(lfc, is.na(sig.label) == FALSE),
             color = "black"
  ) +
  labs( x="Log2FC(Brain/Embryo)", title = "Fly Brain vs Embryo Deseq2") +
  
  ggrepel::geom_text_repel( aes( label = sig.label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 50, alpha = 0.7 ) +
  theme_bw()

ggsave("DEx Brain vs Embryo deseq2.pdf",h=8,w=10)

