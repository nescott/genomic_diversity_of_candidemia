## ---------------------------
## Purpose: Plot intersection of SNPs in series 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------

library(UpSetR)
library(tidyverse)

all_snps <- read.delim("~/tmp4.txt",
                 header = T)

case14 <- all_snps[,c("CHROM", "POS", "MEC131", "MEC132", "MEC133", "MEC134", "MEC324")]

case14 <- case14 %>% 
  right_join(unique_pos)

p <- upset(case14, nsets=8, 
       order.by = "degree", 
       mainbar.y.label = "Shared and Unique SNP Counts", 
       sets.x.label = "Count per isolate",
       main.bar.color = "#88CCEE")

