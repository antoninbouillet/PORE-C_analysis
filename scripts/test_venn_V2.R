library(data.table)
library(dplyr)
library(ggplot2)

rm(list=ls())

baseDir <- "/home/anton/Bureau/PORE-C_repo/data/regions/"

ab <- fread(file = file.path(paste0(baseDir, "all_signif_regions_alpha_STM_500kb.bed")))
astm <- fread(file = file.path(paste0(baseDir, "all_signif_regions_alpha_STM_500kb.bed")))

ab <- ab %>% mutate(ech = ifelse(V6 > 0, "beta", "alpha"))

ggplot(data = ab, aes(x=V2, fill=ech)) +
  geom_bar(alpha=0.5) +
  facet_wrap(~V1)
