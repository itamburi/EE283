library(tidyverse)
library(gridExtra)
library(cowplot)
library(qqman)
library(nycflights13)
library(cowplot)


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



