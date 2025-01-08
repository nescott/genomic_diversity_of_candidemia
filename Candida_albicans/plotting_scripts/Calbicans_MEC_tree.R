## ---------------------------
## Purpose of script: Input, manipulate and plot phylogenetic tree data
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)
## ---------------------------
# Load packages
library(ggtree)
library(paletteer)
library(phangorn)
library(readxl)
library(tidytree)
library(treeio)
library(castor)
library(tidyverse)

species <- "Calbicans"

raxml_file <- "~/umn/data/phylogeny/Calbicans/RAxML_bipartitions.MEC2"

metadata_file <- "~/umn/data/metadata/2024_Calbicans_sorted_patient_info.xlsx"

# Clade colors
color_file <- read.table("ch3/clade_colors.txt")
clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

# Phylogeny file
Candida_raxml <- read.tree(raxml_file) 
Candida_raxml_midpoint <- midpoint(Candida_raxml, node.labels = "support") # keeps branch support

ggtree(Candida_raxml_midpoint) + geom_label2(aes(subset = !isTip, label = label), size =2)

# Clade determination: selecting a branch length (or relative branch length)
longest_dist <- castor::find_farthest_tip_pair(Candida_raxml_midpoint)
Candida_clades <- as_tibble(Candida_raxml_midpoint) %>% 
  mutate(rel_branch = branch.length/longest_dist$distance) %>% # normalize branch length to longest pairwise dist.
  filter(branch.length > 0.03) %>% # set a minimum distance
  filter(!node %in% parent, !grepl("^[MA]", label)) %>% # keep only inner nodes
  rename(clade_node = node) %>% 
  mutate(Clade = c("1", "7", "4", "8A", "8B", "3", "2"))

# Metadata
Candida_metadata <- read_excel(metadata_file)
Candida_metadata[Candida_metadata=="NA"] <- NA
Candida_metadata <- rename(Candida_metadata, label = sample)
Candida_metadata$mec_pt_code <- as.character(Candida_metadata$mec_pt_code)
Candida_metadata$cluster <- as.character(Candida_metadata$pt_cluster)
Candida_metadata$ST <- as.character(Candida_metadata$ST)

# Tree object combining phylo and metadata
Candida_tree <- Candida_raxml_midpoint %>% full_join(Candida_metadata, by = "label")

# Summarize branch lengths by patient
Candida_series_branches <- as_tibble(Candida_tree) %>% 
  filter(!is.na(series)) %>% 
  group_by(series) %>% 
  summarize(max_branch = max(branch.length))

Candida_distances <- castor::get_all_pairwise_distances(Candida_raxml_midpoint)

# Plot the midpoint-rooted tree with patient codes and clade highlighting
# option for labeling a subset with geom_tiplab: aes(subset=(label %in% studies)) aes(subset=(label %in% long_reads))
midpoint_plot <- ggtree(Candida_tree,aes(size = (as.numeric(label) <95 | is.na(as.numeric(label)))))  +
  scale_size_manual(values = c(1, 0.2), guide = "none") +
  geom_tiplab(size = 2, aes(label = mec_isolate_code), align = TRUE, linetype = "dotted", linesize = 0.2) + #
  geom_rootedge() +
  scale_fill_manual(values = clade_colors, guide = "none") +
  geom_treescale(x=0, y=0) +
  geom_cladelab(data = Candida_clades, 
                mapping = aes(node = clade_node, label = Clade),
                align = TRUE,
                offset = 0.0105) +
  geom_hilight(data = Candida_clades, aes(node = clade_node, fill = Clade)) #, align = "both"

ggsave(paste0(Sys.Date(),"_",species,"_midpoint_MEC.pdf"),
       midpoint_plot, bg="white", width = 6, height = 7.5, units = "in")

# Write out sample list in order of tips
write.csv(get_taxa_name(midpoint_plot), 
          "Calbicans_MEC_raxml_midpoint_tips.csv", 
          quote = FALSE,
          row.names = FALSE)