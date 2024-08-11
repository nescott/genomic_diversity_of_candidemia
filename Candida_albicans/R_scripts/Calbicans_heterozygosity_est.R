## ---------------------------
## Purpose: Calculate genome heterozygosity % and plot distribution
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

source("~/umn/thesis/ch3/Calbicans_redcap_summary.R")

# Load packages
library(readxl)
library(tidyverse)
library(writexl)
library(correlation)
library(patchwork)

# Input vars
spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp_depth_paths.txt"
in_patient_data <- "~/umn/data/metadata/2024_Calbicans_sorted_patient_info.xlsx"
ordered_patient_data <-  "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
chr_size <- "~/umn/thesis/ch3/Calbicans_chr_lengths.csv"
save_dir <- "~/umn/images/Calbicans/"
genome_size <- 14324315

# Clade colors
color_file <- read.table("ch3/clade_colors.txt")
clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

# Metadata
pop_data <- read_excel(in_patient_data)
pop_data[pop_data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = TRUE)
pt_order <- pt_order %>%
  left_join(pop_data %>% filter(study=="MEC") %>% select(sample, mec_pt_code, mec_isolate_code))

chrs <- read.csv(chr_size, header = TRUE)

chr_labels <- as_labeller(c("1"="Chr1","2"="Chr2", "3"="Chr3", "4"="Chr4",
                            "5"="Chr5", "6"="Chr6", "7"="Chr7", "8"="ChrR"))

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

snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos)) %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2")) %>% 
  left_join(sample_naming, by = c("sample" = "join_id"))



# Estimate het % for genome and chrs
snp_total <- snp_again %>% 
  group_by(sample) %>%
  summarize(all_snps = sum(snp_count, na.rm = TRUE)) %>%
  mutate(het_percentage = all_snps/genome_size*100)

snp_total <- snp_total %>% 
  left_join(pt_order)%>% 
  left_join(sample_naming, by = c("sample" = "join_id")) %>% 
  mutate(ploidy = case_when(sample == "MEC324" ~ 3,
                            sample == "MEC185" ~ 4,
                            sample %in% c("MEC279", "MEC293", "MEC135", "MEC319", "MEC322") ~ 5,
                            sample %in% c("MEC218", "MEC219", "MEC352", "MEC246", "MEC172", 
                                          "MEC198", "MEC199", "MEC200", "MEC201", "MEC297",
                                          "MEC174") ~ 6,
                            .default = 2))

chr_total <- snp_again %>% 
  group_by(sample, index) %>% 
  summarise(chr_snps = sum(snp_count, na.rm = TRUE)) %>% 
  mutate(chr_het_percent = chr_snps/chrs$length * 100) %>% 
  left_join(pt_order) %>% 
  left_join(sample_naming, by = c("sample" = "join_id"))%>% 
  mutate(ploidy = case_when(sample == "MEC324" ~ 3,
                            sample == "MEC185" ~ 4,
                            sample %in% c("MEC279", "MEC293", "MEC135", "MEC319", "MEC322") ~ 5,
                            sample %in% c("MEC218", "MEC219", "MEC352", "MEC246", "MEC172", 
                                          "MEC198", "MEC199", "MEC200", "MEC201", "MEC297",
                                          "MEC174") ~ 6,
                            .default = 2))

chr_summary <- chr_total %>% group_by(index) %>% 
  summarise(min_het = min(chr_het_percent), max_het = max(chr_het_percent),
            mean_het = mean(chr_het_percent), median_het = median(chr_het_percent))

################################################################################
# Correlation tests

genome_het_gc_corr <- snp_total %>% 
  filter(ploidy==2) %>% 
  left_join(gc, by = "primary_id") %>% 
  select(het_percentage, k, r) %>% 
  correlation()

genome_het_rpmi_ctrl_corr <- snp_total %>% 
  filter(ploidy==2) %>% 
  left_join(control_od, by = "primary_id") %>% 
  select(het_percentage, mean_stationary_k) %>% 
  correlation()

genome_het_smg <- snp_total %>% 
  filter(ploidy==2) %>% 
  left_join(mic_info %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  select(het_percentage, mean_smg) %>% 
  correlation()

chr_het_gc_corr <- chr_total %>% 
  filter(ploidy==2) %>% 
  left_join(gc, by = "primary_id") %>% 
  group_by(index) %>% 
  select(chr_het_percent, k, r) %>% 
  correlation()

chr_het_rpmi_ctrl_corr <- chr_total %>% 
  filter(ploidy==2) %>% 
  left_join(control_od, by = "primary_id") %>% 
  group_by(index) %>% 
  select(chr_het_percent, mean_stationary_k) %>% 
  correlation()

chr_het_smg <- chr_total %>% 
  filter(ploidy==2) %>% 
  left_join(mic_info %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  group_by(index) %>% 
  select(chr_het_percent, mean_smg) %>% 
  correlation()

################################################################################
#Plots

het_all <- ggplot(snp_total, aes(het_percentage)) +
  geom_histogram(color = "grey50", fill = "grey65", binwidth = 0.0105) +
  scale_x_continuous(limits = c(0.0,0.7)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  xlab("Genome heterozygosity, %") +
  ylab("Count of isolates")

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_genome_het.pdf"),
       het_all,
       width = 6,
       height = 4.5,
       units = "in")

het_chr_facet <- chr_total %>% 
  ggplot(aes(chr_het_percent)) +
  geom_histogram(color = "grey50", fill = "grey65", binwidth = 0.01) +
  facet_wrap(~index, nrow = 2, labeller = chr_labels) +
  scale_x_continuous(limits = c(0.0,0.8)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black"),
        ) +
  xlab("Chromosome heterozygosity, %") +
  ylab("Count of isolates")

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_chr_het.pdf"),
       het_chr_facet,
       width = 6,
       height = 4.7,
       units = "in")

gen_chr_het <- (het_all+ theme(plot.margin = margin(b=12)))/
  het_chr_facet + plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_combined_hets.pdf"),
       gen_chr_het,
       width = 6,
       height = 6.2,
       units = "in")
