## ---------------------------
## Purpose of script: Input, manipulate and plot phylogenetic tree data
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)
## ---------------------------
# Load packages
library(tidyverse)
library(paletteer)
library(readxl)
library(ggtree)
library(tidytree)
library(treeio)
library(phangorn)
library(castor)

species <- "Cglabrata"

raxml_file <- "~/umn/data/phylogeny/Cglabrata/RAxML_bipartitions.Cglab_combined"

metadata_file <- metadata_file <- "~/umn/data/metadata/2024_Cglabrata_sorted_patient_info.xlsx"

# Tip colors
patient_colors <- c(paletteer_d("ggsci::default_igv"),
                    "#838b8b",
                    "black",
                    "#cdcdb4",
                    "#155F83")

#st_colors <- c("#59a14f","#ff9da7","#76b7b2","#edc948", "#f28e2b","firebrick","#9c755f", "#4e79a7")
#st_shapes <- c(7,5,6,7,16,17,8,15)
st_shapes <- c(8,12,19)
st_colors <- c("firebrick", "#edc948", "#4e79a7")
st_size <- c(1, 0.3, 0.6,0.7,0.9)

# Phylogeny file
Candida_raxml <- read.tree(raxml_file) 
Candida_raxml_midpoint <- midpoint(Candida_raxml, node.labels = "support") # keeps branch support


# Clade determination: selecting a branch length (or relative branch length)
longest_dist <- castor::find_farthest_tip_pair(Candida_raxml_midpoint)
Candida_clades <- as_tibble(Candida_raxml_midpoint) %>% 
  mutate(rel_branch = branch.length/longest_dist$distance) %>% # normalize branch length to longest pairwise dist.
  filter(rel_branch > 0.05) %>% # set a minimum distance
  filter(!node %in% parent, !grepl("^[[:alpha:]]", label)) %>% # keep only inner nodes
  rename(clade_node = node) 

# Nodes for Gabaldon's clade delineations based on Carrete samples
internal_node <- c(368,385,221,343,324,353,284,256)
gabaldon_clades <- as_tibble(internal_node) %>% 
  mutate(Gclade = c("1", "2", "3", "4", "5", "5", "6", "7"))

# "Clade" labels - sequence types
st_clades <- read_excel(metadata_file, sheet = 2)

# Metadata
Candida_metadata <- read_excel(metadata_file)
Candida_metadata[Candida_metadata=="NA"] <- NA
Candida_metadata <- rename(Candida_metadata, label = sample)
Candida_metadata$mec_pt_code <- as.character(Candida_metadata$mec_pt_code)
Candida_metadata$cluster <- as.character(Candida_metadata$pt_cluster)
Candida_metadata$ST <- as.character(Candida_metadata$ST)

# Tree object combining phylo and metadata
Candida_tree <- Candida_raxml_midpoint %>% full_join(Candida_metadata, by = "label")

# Plot the midpoint-rooted tree with patient codes and clade highlighting , x=0, y=0 , align = TRUE, linetype = "dotted", linesize = 0.2
# option for labeling a subset with geom_tiplab: aes(subset=(label %in% studies)) aes(subset=(label %in% long_reads))
midpoint_plot <- ggtree(Candida_tree,
                        #layout = "circular",
                        aes(size = (as.numeric(label) < 95 | is.na(as.numeric(label)))))  +
  scale_size_manual(values = c(1, 0.3), guide = "none") +
  geom_tiplab(aes(#subset = study %in% c("MEC"), 
                  label = st_label), 
              size=2.2) + #, #align = TRUE, 
              #linetype = "dotted", 
              #linesize = 0.2) + 
  geom_rootedge() +
  geom_treescale(x=0, y=-0.3) +
  geom_tippoint(aes(shape = phylo_label, color = phylo_label), size=0.5) + 
  scale_color_manual(values = st_colors) +
  scale_shape_manual(values = st_shapes) +
  #scale_fill_manual(values = clade_colors, guide = "none") + #name = "Clusters"
  geom_cladelab(data = st_clades, 
                mapping = aes(node = st_nodes, label = st_ID),
                align = TRUE,
                offset = 0.002,
                fontsize = 3
                ) +
  theme(legend.position = c(0.05,0.87), legend.title = element_blank(),
        legend.margin = margin(t =0, b = 0, unit = "cm"),
        legend.text = element_text(margin = margin(l = 0, unit = "cm"))) +
  guides(shape = guide_legend(override.aes = list(size=1.5)))
  #geom_hilight(data = gabaldon_clades, aes(node = value, fill = Gclade)) +
  #scale_fill_manual(values = paletteer_d("ggthemes::Tableau_10"), guide = "none")

ggsave(paste0(Sys.Date(),"_",species,"_midpoint_expanded.pdf"),
       midpoint_plot, bg="white", width = 6, height = 8, units = "in")

# Write out sample list in order of tips
write.csv(get_taxa_name(midpoint_plot), 
          "Cglabrata_expanded_data_raxml_midpoint_tips.csv", 
          quote = FALSE,
          row.names = FALSE)