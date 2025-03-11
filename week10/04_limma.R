here::i_am("04_limma.R")
library(tidyverse)
library(here)
library(edgeR)

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

write.csv(expr2, here("data/batch 3 counts matrix heart only.csv"), row.names = TRUE)






#### 2.0 - Load b2 and b3 counts matrices and combine ####

# batch 2 raw matrix, all tissues, subset for heart
x_b2 = read.csv(here("data/batch 2 counts matrix from old pipeline.csv"), row.names = 1, fileEncoding = "UTF-8-BOM") %>%
  select( matches("Heart"))

# batch 3 raw matrix, heart only
x_b3 = read.csv(here("data/batch 3 counts matrix heart only.csv"), row.names = 1, fileEncoding = "UTF-8-BOM")

names(x_b2) = paste0("b2_", names(x_b2))
names(x_b3) = paste0("b3_", names(x_b3))



#### 3.0 - limma with batches 2 and 3 ####
# https://ucdavis-bioinformatics-training.github.io/2018-September-Bioinformatics-Prerequisites/friday/limma_biomart_vignettes.html

keep = intersect(rownames(x_b2), rownames(x_b3))
expr_h = cbind( x_b2[keep,] , x_b3[keep,] ) 

#### 2.1 - Set up metadata and design matrix for limma ####

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
keep = filterByExpr(dge, design) # note default min count arguments
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
z = topTable(fit2, number=Inf) %>% rownames_to_column(var="gene.id")

z2 = z %>% filter(gene.id %in% c("ISG20","ISG15","PPARA", "CPT2","CRAT", "PPARGC1A", "UCP3","NDUFS1"))




# enrichr
library(enrichR)
dbs1 = c("MSigDB_Hallmark_2020")

gn = z %>% filter(hfhs_vs_nc < 0, adj.P.Val <= 0.05)
enriched1 = enrichr(gn$gene.id, dbs1)
df1 = as.data.frame(enriched1) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
plotEnrich(enriched1[[1]], showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys among downregulated genes relative to NC")

gn = z %>% filter(hfhs_vs_nc > 0, adj.P.Val <= 0.05)
enriched2 = enrichr(gn$gene.id, dbs1)
df2 = as.data.frame(enriched2) %>%
  filter( MSigDB_Hallmark_2020.Adjusted.P.value <= 0.05 )
plotEnrich(enriched2[[1]], showTerms = 30, numChar = 100, y = "Count", orderBy = "P.value", title="Pwys among upregulated genes relative to NC")











