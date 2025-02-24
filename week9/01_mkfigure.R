library(tidyverse)
library(gridExtra)
library(cowplot)
library(qqman)
library(nycflights13)
library(cowplot)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)

#
## Problem 1 - pretty figure for flights dataset ####
P1 <- ggplot(flights, aes(x=distance,y=arr_delay)) +
  geom_point() + 
  geom_smooth() +
  labs(title = "Delays vs Distance", y="Arrival delay (minutes)", x="Flight Distance") +
  theme_bw()

temp_flights <- flights %>%
  group_by(carrier) %>%
  summarize(m_arr_delay = mean(arr_delay,na.rm=TRUE))

P2 <- ggplot(temp_flights, aes(x=carrier,y= m_arr_delay)) +
  geom_bar(stat="identity") +
  labs(title = "Delays by Carrier", y="Mean Arrival delay (minutes)", x="Carrier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

P3 <- ggplot(flights, aes(x=carrier, y=arr_delay)) + 
  geom_boxplot() +
  labs(title = "Distribution of delays by Carrier", y="Arrival delay (minutes)", x="Carrier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

P4 <- ggplot(flights, aes(x=arr_delay)) +
  geom_histogram() +
  labs(title = "Histogram of total delays", y="Count", x="Minutes Delated") +
  theme_bw()

lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,4,4),
             c(3,3,4,4))

tiff("prob1_figure.tiff", width = 7, height = 6, units = "in", res=600)
grid.arrange(P1,P2,P3,P4, layout_matrix = lay)
graphics.off()




## Problem 2 - 4-6 Pane figure from some previous results ####

load("../week7/deseq2_objs.Rdata")

# check dispersions and variance
pp1 = function(){ plotMA( res, ylim = c(-1, 1) ) }
pp2 = function(){ plotDispEsts( dds ) }


# Heatmap of sample distributions
rld = rlog( dds )
mydata = assay(rld)
sampleDists = dist( t( assay(rld) ) )
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pp3 = function(){ heatmap.2( sampleDistMatrix, trace="none", col=colours) }


# heatmap of top discriminating genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
pp4 = function(){
  heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
}


# Volcano plot
lfc = data.frame(res) %>%
  rownames_to_column(var = "gene.id") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(sig.label = ifelse(padj <=0.001, gene.id, NA))

pp5 = ggplot(lfc, aes(log2FoldChange, -log10(padj))) +
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

# plot
title = ggdraw() + draw_label("Dispersion and Variance in DEx Genes \n (Brain over Embryo)", fontface='bold')
z1 = plot_grid(pp1, pp2, labels=c("A","B"))
AB = plot_grid( title, z1, ncol=1, rel_heights = c(.1, 1))

title = ggdraw() + draw_label("Sample Clustering and Top Discriminating Genes", fontface='bold')
z2 = plot_grid(pp3, pp4, labels=c("C","D"))
CD = plot_grid( title, z2, ncol=1, rel_heights = c(.1, 1))

# expand the plot window a lot
row1 = plot_grid(AB, CD)
row2 = plot_grid(NULL, pp5, NULL, nrow = 1, labels =c("", "E", ""), rel_widths = c(.5,1,.5))
figure=plot_grid(row1, row2, ncol=1, rel_heights = c(1,1))
ggsave(plot = figure, "prob2_deseqresults.pdf",h=10, w=20)





## Problem 3 - Make Two Manhattan Plots ####

# mod1
res1 = read.csv("../week8/interaction of treatment to founder at each genome position by two models (neg.log10p).csv") %>%
  mutate(
    SNP = paste0(chr,"_",pos),
    CHR = chr,
    BP = pos,
    P = mod1.signif # notice model 1
  ) %>%
  select(CHR, BP, SNP, P)

chrs = res1 %>%
  distinct(CHR) %>%
  mutate(
    CHR_num = case_when(
      CHR =="chrX" ~ 1,
      CHR =="chr2L" ~ 2,
      CHR =="chr2R" ~ 3,
      CHR =="chr3L" ~ 4,
      CHR =="chr3R" ~ 5
    )
  )

res1_num = res1 %>%
  left_join(.,chrs) %>%
  select(-CHR) %>%
  rename(
    CHR = "CHR_num"
  )


p1 = function(){
  manhattan(res1_num, main = "Significant SNPs from Model1",
          ylim = c(0, 10), cex = 0.6, cex.axis = 0.9,
          col = c("blue4", "orange3"),
          suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08),
          logp=TRUE,
          chrlabs = chrs$CHR)
}


# mod 2
res2 = read.csv("../week8/interaction of treatment to founder at each genome position by two models (neg.log10p).csv") %>%
  mutate(
    SNP = paste0(chr,"_",pos),
    CHR = chr,
    BP = pos,
    P = mod2.signif # notice model 2
  ) %>%
  select(CHR, BP, SNP, P)

chrs = res2 %>%
  distinct(CHR) %>%
  mutate(
    CHR_num = case_when(
      CHR =="chrX" ~ 1,
      CHR =="chr2L" ~ 2,
      CHR =="chr2R" ~ 3,
      CHR =="chr3L" ~ 4,
      CHR =="chr3R" ~ 5
    )
  )

res2_num = res2 %>%
  left_join(.,chrs) %>%
  select(-CHR) %>%
  rename(
    CHR = "CHR_num"
  )


p2 = function(){
  manhattan(res2_num, main = "Significant SNPs from Model 2",
          ylim = c(0, 10), cex = 0.6, cex.axis = 0.9,
          col = c("blue4", "orange3"),
          suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08),
          logp=TRUE,
          chrlabs = chrs$CHR)
}

z = plot_grid(p1, p2, ncol = 1)
print(z)
ggsave("prob3_manhattan plots.pdf",h=9,w=8)

## Problem 4 - Make Manhattan with myManhattan package ####
source("mymanhattan.R")

gg1 = myManhattan(res1_num, chrom.lab =chrs$CHR, graph.title = "Significant SNPs from Model 1")
gg2 = myManhattan(res2_num, chrom.lab =chrs$CHR, graph.title = "Significant SNPs from Model 2")

plot_grid(gg1, gg2, ncol =1)
ggsave("prob4_manhattan from myManhattan.pdf", width=8, height=9)

## Problem 5 - add third pane with comparison of model pvalues ####

merge = read.csv("../week8/interaction of treatment to founder at each genome position by two models (neg.log10p).csv") %>%
  mutate(
    SNP = paste0(chr,"_",pos),
    CHR = chr,
    BP = pos,
    label = ifelse( abs( log10(mod1.signif)-log10(mod2.signif)/mod1.signif)*100 > 50 , SNP,NA)
  )


SNP = ggplot(merge, aes( mod1.signif , mod2.signif) ) +
  labs(
    title = "Comparison of SNP Model Results",
    subtitle = "SNPs with > 50% difference in -log10pval are labeled red",
    x= "-log10 of pvalues from Model 1 ",   y= "-log10 of pvalues from Model 2 "
  ) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = .2, alpha = .4, color = "black") +
  geom_point(size = .2,
             data = subset(merge, is.na(label) == FALSE),
             color = "red"
  ) +
  ggrepel::geom_text_repel( aes( label = label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 100, alpha = 0.7 ) +
  theme_bw()



MP = plot_grid(gg1, gg2, ncol =1)
plot_grid(MP, SNP, nrow=1, rel_widths = c(1.7, 1))
ggsave("prob5_finalfigure.pdf", h=6, w=12)



