
library(tidyverse)
library(cowplot)

# plot the bed files of per bp read coverage
a4=read.delim("/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week4/A4/A4.bedgraph",header = FALSE)
names(a4) =c("chr","start","end","A4")
a5=read.delim("/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week4/A5/A5.bedgraph", header= FALSE)
names(a5) =c("chr","start","end","A5")

z1= merge(a4, a5) %>%
  pivot_longer(cols=c(A4,A5),values_to = "RPKM")
  
gg1 = ggplot(z1,aes(start,RPKM,color=name)) +
  geom_point() +
  geom_line() +
  labs(title="bamCoverage default") +
  scale_color_manual(values=c("blue", "grey80")) +
  geom_vline(xintercept=1904042, color ="red", linetype="dashed") +
  theme_bw()

# plot the bed extended reads
a4=read.delim("/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week4/A4/A4_ext.bedgraph",header = FALSE)
names(a4) =c("chr","start","end","A4")
a5=read.delim("/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week4/A5/A5_ext.bedgraph", header= FALSE)
names(a5) =c("chr","start","end","A5")

z2= merge(a4, a5) %>%
  pivot_longer(cols=c(A4,A5),values_to = "RPKM")

gg2 = ggplot(z2,aes(start,RPKM,color=name)) +
  geom_point() +
  geom_line() +
  labs(title="bamCoverage with --extendReads 500") +
  scale_color_manual(values=c("blue", "grey80")) +
  geom_vline(xintercept=1904042, color ="red", linetype="dashed") +
  theme_bw()

plot_grid(gg1, gg2)
ggsave(h=4,w=16,"/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week4/coveragePlot.png")
