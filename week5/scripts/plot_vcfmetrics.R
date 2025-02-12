
library(ggplot2)
library(cowplot)

data = read.delim("./vcf_filt_metrics.txt", sep = "\t")
data$Step = factor(data$Step, levels =c("Initial VCF", "After Filter 1", "After Filter 2","After Filter 3"))

p1=ggplot(data, aes(Step, SNP_Count)) +
  geom_bar(stat="identity") + 
  labs(x="") +
  theme_bw()

p2=ggplot(data, aes(Step, Mean_DP)) +
  geom_bar(stat="identity") +
  labs(x="") +
  theme_bw()

p3=ggplot(data, aes(Step, Mean_QUAL)) +
  geom_bar(stat="identity") +
  labs(x="") +
  theme_bw()

p4=ggplot(data, aes(Step, TiTv_Ratio)) +
  geom_bar(stat="identity") +
  labs(x="") +
  theme_bw()

plot_grid(p1,p2,p3,p4, nrow=2)
ggsave("snps after 3 filtering steps.png",h=6,w=8)
