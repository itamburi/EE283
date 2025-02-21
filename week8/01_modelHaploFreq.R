library(tidyverse)
library(broom)
getwd()

## 0.0 - Handy tricks from class notes ####

# (a) approximate joins:
closest_gene = hits %>% left_join(genes, join_by(chr, closest(pos)))
#This command is saying to left_join using the shared columns chr and pos (only), with an exact match for chr, but only the closest pos in the genes table for each row in hits.


# (b) apply a function to grouped portions of the data with nest, map and broom::glance

mySimpleModel <- function(df){
  lm(freq ~ treat + rep %in% treat, data=df)
}

library(broom)
myresult = malathion %>%
  group_by(chr, pos) %>%
  nest() %>%
  # Map applies the function to the "data" column of nest, which contains a df in each cell
  mutate(myMod = map(data, mySimpleModel)) %>%
  # broom::glance extracts sumamry statistincs from myMod
  mutate(myGlance = map(myMod, broom::glance)) %>%
  unnest(myGlance, .drop = TRUE)

# Here the map function returns an “object” that results from a linear model
# fit to the data.  The broom function glance extracts several useful summary
# statistics from the linear model, and then saves them as this object called
# myGlance.  Like “data” and “myMod”, “myGlance” is a funny object in the tibble
# we cannot see (without extracting and looking at individual instances).  In
# the final step we use the “unnest” function to extract these summaries as new
# variables (that are visible) in the resulting tibble. The utility of broom
# seems to be that the objects created by it (for example via glance) are
# intrinsically un-nestable. Whereas if you just tried to unnest myMod it is
# not clear that action has any meaning.


## 1.0 - Homework Problems ####

mal = read_tsv("allhaps.malathion.200kb.txt.gz")

# select the first position in the genome
mal2 = mal %>% filter(chr=="chrX" & pos==316075)
mal2
levels(as.factor(mal2$pool))
levels(as.factor(mal2$founder))
# 32 obs per position in genome
# The frequencies of 8 haplotypes (founders) are measured in four populations: mcF, mcM, msF, msM; male or female
# obtained from a cage prior to (=“c”ontrol) or after “s”election with malathion (a pesticide)
# The idea of the experiment is that protective haplotypes at different locations in the genome should show enrichment in the malathion treated flies.

# fitting a model for change at a position in the genome:
mal2 = mal2 %>% mutate(treat=str_sub(pool,2,2))
anova(lm( asin(sqrt(freq)) ~ treat + founder + treat:founder, data=mal2))
# The p-value of interest is the two-way interaction between treat and founder
# (which is the same as the effect of treatment within founder after the effect of treatment is removed) = 0.09027.



# (1) - fit this model and extract the -log10(p) at every location in the genome using group_by, nest, and map.
mal = read_tsv("allhaps.malathion.200kb.txt.gz") %>%
  mutate(treat=str_sub(pool,2,2))
names(mal)

mymod1 <- function(df){
  anova(lm( asin(sqrt(freq)) ~ treat + founder + treat:founder, data=df))
}

extr_pval1 <- function(df) {
  pval <- df$`Pr(>F)` # Extracts p-values for all terms of model
  return(pval[3])   # treat:founder is in row 3
}

aov = mal %>%
  filter(
    is.na(freq) == FALSE
  ) %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(res = map(data, mymod1)) %>%
  mutate(
    p = map(res, extr_pval1),
    neg.log10p = -log10(as.numeric(p))
         ) %>%
  select(!c("data", "res")) %>%
  unnest()

write.csv(aov, "model1 interaction of treatment to founder at each genome position.csv")


# (2) - Repeat scan with alternative model.

mymod2 <- function(df){
  anova(lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data=df))
}

extr_pval2 <- function(df) {
  pval <- df$`Pr(>F)` # Extracts p-values for all terms of model
  return(pval[2])   # treat:founder is in row 2
}


aov = mal %>%
  filter(
    is.na(freq) == FALSE
  ) %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(res = map(data, mymod2)) %>%
  mutate(
    p = map(res, extr_pval2),
    neg.log10p = -log10(as.numeric(p))
  ) %>%
  select(!c("data", "res")) %>%
  unnest()
write.csv(aov, "model2 interaction of treatment to founder at each genome position.csv")


# (3) - merge the two scans into a single table

mod1 = read.csv("model1 interaction of treatment to founder at each genome position.csv") %>%
  select(chr, pos, neg.log10p) %>%
  rename(mod1.signif = "neg.log10p")

mod2 = read.csv("model2 interaction of treatment to founder at each genome position.csv") %>%
  select(chr, pos, neg.log10p) %>%
  rename(mod2.signif = "neg.log10p")

idx1 = mod1 %>%
  distinct(chr, pos) %>%
  mutate( idx = paste0(chr,pos) )

idx2 = mod2 %>%
  distinct(chr, pos) %>%
  mutate( idx = paste0(chr,pos) )

identical(idx1$idx, idx2$idx)

merge = mod1 %>% left_join(., mod2)

write.csv(merge, "interaction of treatment to founder at each genome position by two models (neg.log10p).csv", row.names = FALSE)




