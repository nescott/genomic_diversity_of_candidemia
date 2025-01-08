## ---------------------------
## Purpose: Combine het SNP counts and generate heatmap across chrs for a group of isolates
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(writexl)

spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp_depth_paths.txt"
ordered_patient_data <- "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
in_patient_data <- "~/umn/data/metadata/2024_Calbicans_sorted_patient_info.xlsx"
save_dir <- "~/umn/images/Calbicans/"

binned_color_scale <- c("white", paletteer_c("grDevices::Turku", 30))

# Read in SNP counts for all samples
snp_files <- scan(spreadsheet_list, what=character())

genome_snp <- read_xlsx(snp_files[1]) %>%
  select(index, pos, snp_count)
names(genome_snp)[names(genome_snp)=="snp_count"] <-str_extract(snp_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:100){
  new_snp <- read_xlsx(snp_files[i]) %>%
    select(index, pos, snp_count)
  names(new_snp)[names(new_snp)=="snp_count"] <-str_extract(snp_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_snp <- genome_snp %>%
    left_join(new_snp, by = join_by(index,pos))
}

# For plotting: regular intervals for x-axis
genome_snp$x_num <- seq.int(nrow(genome_snp))

# Finish tidying data
snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos, x_num))

snp_again <- snp_again %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2"))

# For plotting: chr outlines and tick marks
chrs <- genome_snp %>%
  group_by(index) %>%
  summarise(border_start=min(x_num),
            border_stop=max(x_num),
            tick=min(x_num) + (max(x_num)-min(x_num))/2
            )

# Metadata, sort samples by tree order
pop.data <- read_excel(in_patient_data)
pop.data[pop.data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = TRUE)
pt_order <- pt_order %>%
  left_join(pop.data %>% filter(study=="MEC") %>% select(sample, mec_pt_code, mec_isolate_code))

################################################################################
# Plot
p <- snp_again %>%
  mutate(clustered_samples = fct_relevel(sample, rev(pt_order$sample))) %>%
  ggplot(aes(x=x_num, y=clustered_samples)) +
  geom_tile(aes(fill=snp_count, colour = snp_count), linewidth=0.27) +
  scale_color_stepsn(na.value = "white", 
                     limits = c(0,300),
                     n.breaks = 31,
                     colors = binned_color_scale,
                     guide = "none") +
  scale_fill_stepsn(na.value = "white", 
                    limits = c(0,300),
                    n.breaks = 31,
                    labels = c("0", rep("",9), "100", rep("", 9), "200", rep("", 9), "300"),
                    colors = binned_color_scale,
                    name = "Heterozygoups SNPs per 5kb") +
  scale_y_discrete(labels = rev(pt_order$mec_isolate_code)) +
  theme_minimal() +
  geom_rect(data=chrs, aes(group=index, xmin=border_start, xmax=border_stop, ymin=0.5, ymax=Inf),
                            fill=NA, inherit.aes=FALSE,colour = "grey26", linejoin = "round") +
  scale_x_continuous(expand = c(0,0),
                      breaks = chrs$tick,
                      labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "ChrR")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=5.4, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.position = "bottom", 
        legend.title.position = "left",
        legend.title = element_text(size = 8, color = "black"),
        plot.margin = margin(0,10.5,5.5,5.5,"pt"),
        legend.text = element_text(angle = 90, size = 6, color = "black")) +
  annotate("segment", x=0, xend = 2860, y= 19.5, yend = 19.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 20.5, yend = 20.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 34.5, yend = 34.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 35.5, yend = 35.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 36.5, yend = 36.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 37.5, yend = 37.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 43.5, yend = 43.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 46.5, yend = 46.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 54.5, yend = 54.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 19.5, yend = 19.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 59.5, yend = 59.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 63.5, yend = 63.5, linewidth = 0.2, colour = "grey39") 

# Troubleshooting: what are the axis positions per sample?
#details <- ggplot_build(p)

ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_LOH_heatmap.pdf"), 
       p, 
       #device=png, 
       #dpi=300, 
       #bg="white",
       width = 6, 
       height = 6, 
       units="in")

