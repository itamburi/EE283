
library(ggplot2)
library(cowplot)
library(pheatmap)


# Prob #3 - plot all of the snps across genotypes
gen = read.delim("X_1Mb_snps.012", sep = "\t", header=FALSE)
gen = as.data.frame(t(gen[,-1]))
pos = read.delim("X_1Mb_snps.012.pos", sep = "\t", header=FALSE)

rownames(gen) = pos$V2
names(gen) = c("geneotype1","geneotype2","geneotype3","geneotype4")

color_palette <- c("green", "yellow", "red")  
breaks <- c(-0.5, 0.5, 1.5, 2.5)

pheatmap(
    gen,
    main="Biallelic distribution across first 1MB of Chrx\n Green = ref/ref\n Yellow = ref/alt\n Red = alt/alt",
    cluster_rows = FALSE, cluster_cols = FALSE,
    fontsize_row = 4,
    color = color_palette, breaks = breaks
    )

ggsave("genotype_heatmap.png", plot = grid::grid.grab(), width = 8, height = 12)


# Prob #4 - only strains where two are 0/0 and other two are 1/1; snp at frequency of 50%
library(tidyverse)

# this is trivial to perform in R. But It can also be done with awk; see scripts/08_awkhomozyg.sh
filt = gen %>%
    filter(!if_any(everything(), ~ . == 1)) %>%
    mutate(sum = geneotype1+geneotype2+geneotype3+geneotype4) %>%
    filter(sum==4) %>%
    mutate(sum = NULL)

pheatmap(
    filt,
    main="SNPs with 50% freq isogenicity in chrX\n Green = ref/ref\n Yellow = ref/alt\n Red = alt/alt",
    cluster_rows = FALSE, cluster_cols = FALSE,
    fontsize_row = 4,
    color = color_palette, breaks = breaks
)

ggsave("genotype_heatmap_freq50.png", plot = grid::grid.grab(), width = 6, height = 6)



