## ---------------------------
## Purpose: Calculate and plot pairwise differing sites or alleles
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
#options(scipen = 999)
## ---------------------------
## load packages
library(vcfR)
library(poppr)
library(tidyverse)
library(readxl)
library(rstatix)
library(scales)
library(patchwork)

source("redcap_MIC_summary.R")

species <- "C. albicans"
in_variant_file <- "~/umn/data/variants/Calbicans/Calbicans_MEC_bwa_filtered_annotated.vcf.gz"
in_patient_data <- "~/umn/data/metadata/2024_Calbicans_sorted_patient_info.xlsx"
ordered_patient_data <-  "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
save_dir <- "~/umn/images/Calbicans/"

species_ploidy <- 2

# Clade colors
color_file <- read.table("ch3/clade_colors.txt")
clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

## create genlight object from vcf file (can be .gz) and separate sample population file
vcf <- read.vcfR(in_variant_file)

pop.data <- read_excel(in_patient_data)
pop.data[pop.data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = TRUE)
pt_order <- pt_order %>%
  left_join(pop.data %>% filter(study=="MEC") %>% select(sample, mec_pt_code, mec_isolate_code))

all(colnames(vcf@gt)[-1] %in% pop.data$sample[1:100])

gl_species <- vcfR2genlight(vcf)
ploidy(gl_species) <- species_ploidy

# Get distance matrix
# setting differences_only to TRUE matches plink 1.9 --genome full output (IBS0 + IBS1 per row)
# differences_only FALSE counts differences in alleles, not just genotypes
intra_species_false <- bitwise.dist(gl_species, percent = FALSE, mat = TRUE, differences_only = FALSE)
intra_species_true <- bitwise.dist(gl_species, percent = FALSE, mat = TRUE, differences_only = TRUE)

# Pivot longer, add names
intra_species_false <- intra_species_false[pt_order$sample, pt_order$sample]

snp_dist <- intra_species_false %>%
  pull_upper_triangle(diagonal = TRUE) %>% 
  as_tibble() %>%
  rename(Var1 = rowname) %>% 
  mutate(Var1 = rownames(intra_species_false))
 
snp_dist <- snp_dist %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>% 
  filter(value !="") %>% 
  mutate(value = as.integer(value))

snp_dist <- snp_dist %>% 
  left_join(pt_order %>% select(sample, mec_pt_code, cluster), by=c("Var1" = "sample")) %>%
  rename(pt_code_1 = mec_pt_code, cluster1 = cluster)

# Get cluster info
snp_dist <- snp_dist %>% 
  left_join(pt_order %>% select(sample, mec_pt_code, cluster), by=c("Var2" = "sample")) %>%
  rename(pt_code_2 = mec_pt_code, cluster2 = cluster)

# Get patient info for adding days
sample_info <- sample_info %>% 
  mutate(primary_id = case_when(!is.na(secondary_id) ~ secondary_id, .default = primary_id))

# Within-strain vals and summaries
intra_strain_dist <- snp_dist %>% 
  filter(pt_code_1==pt_code_2, Var1!=Var2) %>% 
  distinct(value,pt_code_1,cluster1,pt_code_2,cluster2, .keep_all = TRUE) %>% 
  left_join(sample_info %>% select(primary_id, relative_days), by=c("Var1" = "primary_id")) %>% 
  rename(relative_days1 = relative_days) %>% 
  left_join(sample_info %>% select(primary_id, relative_days), by=c("Var2" = "primary_id")) %>% 
  rename(relative_days2 = relative_days)
  
intra_strain_dist <- intra_strain_dist %>% 
  mutate(days_apart = abs(relative_days1 - relative_days2)) %>% 
  arrange(pt_code_1, value, days_apart)

serial_intra_dist <- intra_strain_dist %>% 
  filter(days_apart < 20)

# Isolate collection ranges mean don't bother reporting this
correlation::correlation(serial_intra_dist, select = c("value", "days_apart"))
ggplot(intra_strain_dist, aes(x = value, y = days_apart, color = cluster1)) + 
  geom_point() #+
  #scale_color_manual(values = paletteer_d("ggthemes::Tableau_20"))


clade_pt_intra_dist_summary <- intra_strain_dist %>% 
  filter(!cluster1 %in% c("S2", "S3", "S4", "S5")) %>% 
  group_by(cluster1, pt_code_1) %>% 
  summarise(min_dist = min(value), max_dist = max(value)) %>% 
  arrange(cluster1, max_dist)

cladewise_intra_dist_summary <- clade_pt_intra_dist_summary %>% 
  filter(cluster1 != "S1") %>% 
  group_by(cluster1) %>% 
  summarise(min_intra_dist = min(min_dist), max_intra_dist = max(max_dist), pt_count = n()) %>% 
  arrange(cluster1) 

# Between strain vals and summary
inter_strain_dist <- snp_dist %>% 
  filter(pt_code_1!=pt_code_2) %>% 
  arrange(cluster1,value)

cladewise_inter_dist_summary <- inter_strain_dist %>% 
  filter(!cluster1 %in% c("S1", "S2", "S3", "S4", "S5")) %>% 
  group_by(cluster1,cluster2) %>%
  summarise(min_inter_dist = min(value), max_inter_dist = max(value)) %>% 
  filter(cluster1==cluster2)

cladewise_all_dist_summary <- cladewise_intra_dist_summary %>% 
  left_join(cladewise_inter_dist_summary) %>% 
  mutate(ratio_min = round(min_inter_dist/max_intra_dist, digits = 2), 
         ratio_max = round(max_inter_dist/max_intra_dist, digits = 2))

# Save dist counts
write.csv(snp_dist, "~/umn/data/variants/Calbicans/2023_Calbicans_MEC_pairwise_snp_count.csv", row.names = FALSE)
#write.csv(intra_strain_dist, "~/umn/data/variants/Calbicans/2023_Calbicans_MEC_intra_strain_snps.csv", row.names = FALSE)
#write.csv(inter_strain_dist, "~/umn/data/variants/Calbicans/2023_Calbicans_MEC_inter_strain_snps.csv", row.names = FALSE)
write.csv(cladewise_all_dist_summary, "~/umn/data/variants/Calbicans/2023_Calbicans_cladewise_all_dists.csv", row.names = FALSE)

# Plot heatmap of pairwise SNP distances
snp_matrix <- ggplot(snp_dist, aes(Var1, Var2)) +
  geom_tile(aes(fill = value, colour = value),linewidth = 0.27) +
  scale_color_gradient(na.value = "white", low = "white", high = "#B22222",guide = "none") +
  scale_fill_gradient(na.value = "white", low = "white", high = "#B22222", name = "Pairwise SNP distance",
                      labels = label_scientific()) +
  scale_x_discrete(limits = pt_order$sample, labels = pt_order$mec_isolate_code, position = "top") +
  scale_y_discrete(limits = pt_order$sample, labels = pt_order$mec_isolate_code) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.text.y = element_text(size = 5, color = "black")) +
  theme(axis.text.x= element_text(angle = 90, size = 5.5, color = "black")) +
  theme(legend.position = "inside", 
        legend.direction = "horizontal",
        #legend.title.position = "top",
        legend.justification = c("center", "bottom"),
        legend.text = element_text(angle = 90, size = 5.5))

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_pairwise_SNP_heatmap.pdf"),
       snp_matrix,
       bg = "white",
       width = 6,
       height = 6,
       units = "in")

# Plot cladewise SNP distributions
intra_range <- ggplot(snp_dist %>% filter(pt_code_1 == pt_code_2, !cluster1 %in% c("S2", "S3", "S4", "S5")), 
                      aes(x = cluster1, y=value, fill = cluster1)) + 
  geom_violin(alpha = 0.6) + 
  scale_fill_manual(values=clade_colors, guide = "none") + 
  theme_bw() + 
  #scale_y_continuous(limits = c(0,10000)) +
  ylab("Intra-strain SNP distance") + 
  xlab("Phylogenetic cluster")

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_intra_strain_snp_distribution.pdf"),
       intra_range,
       bg = "white",
       width = 6,
       height = 4.5,
       units = "in")

inter_range <- ggplot(snp_dist %>% filter(pt_code_1 != pt_code_2), 
                      aes(x = cluster1, y=value, fill = cluster1)) + 
  geom_violin(alpha = 0.6) + 
  scale_fill_manual(values=clade_colors, guide = "none") + 
  theme_bw() + 
  ylab("Inter-strain SNP distance") + 
  xlab("Phylogenetic cluster")

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_inter_strain_snp_distribution.pdf"),
       inter_range,
       bg = "white",
       width = 6,
       height = 4.5,
       units = "in")
