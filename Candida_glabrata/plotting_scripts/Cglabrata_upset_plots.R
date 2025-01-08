## ---------------------------
## Purpose: Plot intersection of SNPs in series 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------

library(UpSetR)
library(tidyverse)

snps <- read.csv("~/umn/data/variants/Cglabrata/Cglabrata_MEC_BB_SNPS.csv",
                 header = T)

names(snps)[names(snps)=="MEC328"] <- "54-1"
names(snps)[names(snps)=="MEC329"] <- "54-2"
names(snps)[names(snps)=="MEC330"] <- "54-3"
names(snps)[names(snps)=="MEC331"] <- "54-4"
names(snps)[names(snps)=="MEC332"] <- "54-5"
names(snps)[names(snps)=="MEC333"] <- "54-6"
names(snps)[names(snps)=="MEC334"] <- "54-7"
names(snps)[names(snps)=="MEC335"] <- "54-8"

res_colors <- c(rep("#004488", 5), "#bb5566", rep("#004488", 2))

p <- upset(snps, nsets=8, 
       order.by = "degree", 
       mainbar.y.label = "Shared and Unique SNPs", 
       sets.x.label = "Within-series SNP count",
       sets.bar.color = res_colors,
       keep.order = T, 
       sets = rev(c("54-1", "54-2", "54-3", "54-4", "54-5", "54-6", "54-7", "54-8")),
       text.scale = c(1.5,1.2,1.5,1,1.2,1.2),
       queries = list(list(query = intersects, params = list("54-3"), active = T,
                           color = "#bb5566")))


pdf("Cglabrata_case54_upset.pdf", width = 6, height = 5.5, bg = "white")
print(p)
dev.off()
