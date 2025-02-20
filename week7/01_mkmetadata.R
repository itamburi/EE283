library(tidyverse)

getwd()
mytab = read_tsv("../datalinks/RNAseq/RNAseq384_SampleCoding.txt")
mytab

mytab2 <- mytab %>%
  select(RILcode, TissueCode, Replicate, FullSampleName)	
table(mytab2$RILcode)
table(mytab2$TissueCode)
table(mytab2$Replicate)


mytab2 <- mytab %>%
  select(RILcode, TissueCode, Replicate, FullSampleName) %>%
  filter(RILcode %in% c(21148, 21286, 22162, 21297, 21029, 22052, 22031, 21293, 22378, 22390)) %>%
  filter(TissueCode %in% c("B", "E")) %>%
  filter(Replicate == "0")

for(i in 1:nrow(mytab2)){
  cat("/dfs6/pub/itamburi/ee283/week3/prob5/RNAbam/",mytab2$FullSampleName[i],".sort.bam\n",file="shortRNAseq.names.txt",append=TRUE,sep='')
}
write_tsv(mytab2,"shortRNAseq.txt")
