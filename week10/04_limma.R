here::i_am("04_limma.R")
library(tidyverse)
library(here)
library(edgeR)
library(enrichR)
library(cowplot)

# *** Batch 3 is the most recent RNA-seq from LDLR-/- animals fed a western diet for 3 or 5 months. These are the reads we aligned in this project
# *** Batch 2 was produced previously on Normal genetic animals fed a regular chow diet or a western-diet
# *** All together, 4 cohorts: Chow, HFHS, LDLR-3mo, LDLR-5mo
# *** For this analysis we will compare heart tissue gene expression between the 4 cohorts using Limma
#
#### 1.0 - Assemble counts from batch 3 pipeline into a matrix ####

# list files and merge 
cfiles = list.files(here("counts"))
counts_list = list()

for (i in 1:length(cfiles)) {
  
  path = paste0("counts/", cfiles[i])
  c = read.delim(here(path), header=FALSE)
  names(c) = c("gene.id", cfiles[i])
  rownames(c) = c$gene.id
  c$gene.id = NULL
  counts_list[[i]] = c
  
}

# make a matrix
matrix = bind_cols(counts_list)
names(matrix) = gsub("_counts.txt","", names(matrix))
matrix$gene.id = row.names(matrix)

# merge with the metadata for sample IDs
matrix_long = reshape2::melt(matrix)

metadata = read.csv(here("data/batch 3 metadata.csv"))

expr = left_join(matrix_long, metadata, by = c("variable" = "novo_id"))

expr2 = expr %>%
  select(internal_id, gene.id, value) %>%
  pivot_wider(names_from = internal_id, values_from = value) %>%
  column_to_rownames(var = "gene.id")

#write.csv(expr2, here("data/batch 3 counts matrix.csv"), row.names = TRUE)
write.csv(expr2, here("data/batch 3 counts matrix heart only.csv"), row.names = TRUE)






#### 2.0 - Load b2 and b3 counts matrices and combine ####

# batch 2 raw matrix, all tissues, subset for heart
x_b2 = read.csv(here("data/batch 2 counts matrix from old pipeline.csv"), row.names = 1, fileEncoding = "UTF-8-BOM") %>%
  select( matches("Heart"))

# batch 3 raw matrix, heart only
x_b3 = read.csv(here("data/batch 3 counts matrix.csv"), row.names = 1, fileEncoding = "UTF-8-BOM")

names(x_b2) = paste0("b2_", names(x_b2))
names(x_b3) = paste0("b3_", names(x_b3))



#### 3.0 - limma with batches 2 and 3 ####
# https://ucdavis-bioinformatics-training.github.io/2018-September-Bioinformatics-Prerequisites/friday/limma_biomart_vignettes.html

keep = intersect(rownames(x_b2), rownames(x_b3))
expr_h = cbind( x_b2[keep,] , x_b3[keep,] ) %>%
  select(matches("_Heart"))

#### 4.0 - Set up metadata and design matrix for limma ####

metadata = data.frame(sample = names(expr_h) ) %>%
  separate(sample, into=c("batch","animal","tissue"), sep="_", remove = FALSE) %>%
  mutate(
    cohort = case_when(
      grepl("nc", animal) ~ "chow",
      grepl("hf", animal) ~ "hfhs",
      grepl("ldl0[0-9]|ldl10", animal) ~ "3mo", # distinguish 3 and 5 month for anova
      grepl("ldl[11-20]", animal) ~ "5mo"
    ),
    cohort = factor(cohort, levels =c("chow","hfhs","3mo","5mo"))
  ) %>%
  column_to_rownames(var="sample")

design <- model.matrix(~cohort+batch, data=metadata)


# DGEList object is a container for counts, normalization factors, and library sizes.
dge = DGEList(counts=expr_h)

#remove rows that consistently have zero or very low counts
keep = filterByExpr(dge, design, min.count = 20, min.total.count = 30) # note default min count arguments
dge = dge[keep,,keep.lib.sizes=FALSE]

filtnames = rownames(dge$counts)


# It is usual to apply scale normalization to RNA-seq read counts,
# and the TMM normalization method in particular has been found to perform well in comparative studies
dge = calcNormFactors(dge)
summary(dge$samples$lib.size)

# When the library sizes are quite variable between samples, then the voom approach is theoretically more powerful than limma-trend
v = voom(dge, design, plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)
z = topTable(fit) %>% rownames_to_column(var="gene.id")

contrasts <- makeContrasts(
  hfhs_vs_nc = cohorthfhs - Intercept,
  cohort3mo_vs_nc = cohort3mo - Intercept,
  cohort5mo_vs_nc = cohort5mo - Intercept,
  levels = design
)


fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
top = topTable(fit2, number=Inf) %>% rownames_to_column(var="gene.id")

gn_up = top %>%
  filter(
    hfhs_vs_nc > 0 & hfhs_vs_nc < cohort3mo_vs_nc & cohort3mo_vs_nc < cohort5mo_vs_nc,
    #hfhs_vs_nc > 0 & hfhs_vs_nc < cohort3mo_vs_nc & hfhs_vs_nc < cohort5mo_vs_nc,
    #cohort3mo_vs_nc > 0 & cohort3mo_vs_nc < cohort5mo_vs_nc,
    adj.P.Val <= 0.05
    )
gn_dn = top %>%
  filter(
    hfhs_vs_nc < 0 & hfhs_vs_nc > cohort3mo_vs_nc & cohort3mo_vs_nc > cohort5mo_vs_nc,
    adj.P.Val <= 0.05
  )


#### 5.0 - pathway enrichment analysis ####
dbs1 = c("MSigDB_Hallmark_2020")

# HFHS vs control
hfhs_up = top %>% filter(hfhs_vs_nc > 0, adj.P.Val <= 0.05)
enriched1 = enrichr(hfhs_up$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee1 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys up in HFHS vs NC")

hfhs_dn = top %>% filter(hfhs_vs_nc < 0, adj.P.Val <= 0.05)
enriched1 = enrichr(hfhs_dn$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee2 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys down in HFHS vs NC")

## LDLR-3mo vs control
ldl3mo_up = top %>% filter(cohort3mo_vs_nc > 0, adj.P.Val <= 0.05)
enriched1 = enrichr(ldl3mo_up$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee3 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys up in 3mo-LDLR-/- vs NC")

ldl3mo_dn = top %>% filter(cohort3mo_vs_nc < 0, adj.P.Val <= 0.05)
enriched1 = enrichr(ldl3mo_dn$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee4 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys down in 3mo-LDLR-/- vs NC")


## LDLR-5mo vs control
ldl5mo_up = top %>% filter(cohort5mo_vs_nc > 0, adj.P.Val <= 0.05)
enriched1 = enrichr(ldl3mo_up$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee5 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys up in 5mo-LDLR-/- vs NC")

ldl5mo_dn = top %>% filter(cohort5mo_vs_nc < 0, adj.P.Val <= 0.05)
enriched1 = enrichr(ldl3mo_dn$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
ee6 = plotEnrich(enriched1[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys down in 5mo-LDLR-/- vs NC")




#### 6.0 - Volcano plots ####

lfc1 = top %>%
  select("gene.id", "adj.P.Val", "hfhs_vs_nc") %>%
  rename(lfc = "hfhs_vs_nc") %>%
  mutate(sig.label = ifelse(adj.P.Val <= 0.001 & abs(lfc) >1.5, gene.id, NA))

vol1 = ggplot(lfc1, aes(lfc, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", alpha = .5) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", alpha = .5) +
  geom_point(size = .2, alpha = .4, color = "grey80") +
  geom_point(size = .2,
             data = subset(lfc1, is.na(sig.label) == FALSE),
             color = "black"
  ) +
  labs( x="Log2FC(HFHS/NC)", title = "HFHS vs NC DEx in Heart") +
  
  ggrepel::geom_text_repel( aes( label = sig.label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 50, alpha = 0.7 ) +
  theme_bw()



lfc2 = top %>%
  select("gene.id", "adj.P.Val", "cohort3mo_vs_nc") %>%
  rename(lfc = "cohort3mo_vs_nc") %>%
  mutate(sig.label = ifelse(adj.P.Val <= 0.001 & abs(lfc) >1.5, gene.id, NA))

vol2 = ggplot(lfc2, aes(lfc, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", alpha = .5) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", alpha = .5) +
  geom_point(size = .2, alpha = .4, color = "grey80") +
  geom_point(size = .2,
             data = subset(lfc2, is.na(sig.label) == FALSE),
             color = "black"
  ) +
  labs( x="Log2FC(3mo LDLR-KO/NC)", title = "3mo LDLR-KO vs NC DEx in Heart") +
  
  ggrepel::geom_text_repel( aes( label = sig.label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 50, alpha = 0.7 ) +
  theme_bw()


lfc3 = top %>%
  select("gene.id", "adj.P.Val", "cohort5mo_vs_nc") %>%
  rename(lfc = "cohort5mo_vs_nc") %>%
  mutate(sig.label = ifelse(adj.P.Val <= 0.001 & abs(lfc) >1.5, gene.id, NA))

vol3 = ggplot(lfc3, aes(lfc, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", alpha = .5) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", alpha = .5) +
  geom_point(size = .2, alpha = .4, color = "grey80") +
  geom_point(size = .2,
             data = subset(lfc2, is.na(sig.label) == FALSE),
             color = "black"
  ) +
  labs( x="Log2FC(5mo LDLR-KO/NC)", title = "5mo LDLR-KO vs NC DEx in Heart") +
  
  ggrepel::geom_text_repel( aes( label = sig.label ),
                            vjust = 1.0,
                            box.padding = 0.5,
                            size = 2.0,
                            max.overlaps = 50, alpha = 0.7 ) +
  theme_bw()


#### 7.0 - Make the figure ####
a = plot_grid(ee2, vol1, ee1, rel_widths = c(1,1.5,1), nrow = 1)
b = plot_grid(ee4, vol2, ee3, rel_widths = c(1,1.5,1), nrow = 1)
c = plot_grid(ee6, vol3, ee5, rel_widths = c(1,1.5,1), nrow = 1)

f = plot_grid(a,b,c, ncol=1, labels = c("A","B","C"))
ggsave(plot=f, "Differentially expresse pathways across cohorts.pdf",h=12,w=20)

