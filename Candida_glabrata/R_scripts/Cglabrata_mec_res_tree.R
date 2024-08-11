source("~/umn/thesis/ch2/redcap_MIC_summary.R", echo=TRUE)
source("~/umn/thesis/ch4/Cglabrata_MEC_tree.R")

library(patchwork)
resistant_isolates <- resistant_isolates %>% rename(label = primary_id)
resistant_isolates <- resistant_isolates %>% filter(genus_species=="C. glabrata")
resistant_isolates <- resistant_isolates %>% 
  mutate(res_cat = case_when(eucast_breakpoint_micafungin=="R" & is.na(eucast_breakpoint_fluconazole) ~ "mcf_r",
                             eucast_breakpoint_fluconazole == "R" & is.na(eucast_breakpoint_micafungin) ~ "flc_r",
                             eucast_breakpoint_micafungin== "R"& eucast_breakpoint_fluconazole == "R" ~ "mdr"))

# more metadata
Candida_metadata <-Candida_metadata %>% left_join(resistant_isolates %>% select(label, res_cat)) %>% 
  filter(study=="MEC") %>% 
  mutate(res_cat = case_when(is.na(res_cat) ~ "s",
                             .default = res_cat))

Candida_metadata$res_cat[Candida_metadata$label=="MEC014"] <- "s"

sample_order <- read_csv("~/umn/data/metadata/Cglabrata_MEC_raxml_midpoint_tips.csv") %>% 
  rename(label = sample) %>% 
  left_join(Candida_metadata %>% select(label, mec_isolate_code, ST))

ch4_table <- read_xlsx("~/Google Drive/My Drive/dissertation/04 Genomic plasticity of Candida glabrata bloodstream isolates/Ch4_supp_data.xlsx",
                       sheet = 2)
ch4_table <- ch4_table %>%
  rename(label = sample) %>% 
  mutate(cnv = case_when(`Mean relative depth` < 0.5 ~ "del",
                         `Mean relative depth` > 1.4 ~ "amp")) %>% 
  left_join(Candida_metadata %>% select(label, st))

try1 <- ch4_table %>% 
  select(label, Locus, cnv) %>%
  group_by(label) %>%
  pivot_wider(names_from = Locus, values_from = cnv, values_fill = NA)

cnvs <- ggplot(ch4_table, 
               aes(x = Locus, y = label, fill = cnv, color = "lightgrey")) + 
  geom_raster(hjust = 0, vjust =0)  +
  scale_y_discrete(limits = rev(sample_order$label), 
                   labels = rev(sample_order$mec_isolate_code), 
                   #position = "right", 
                   expand = expansion(add = c(0,0))) +
  theme(axis.title =element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, color = "black", size = 4, vjust = -0.8),
        axis.text.y = element_text(color = "black", size = 5.5, vjust = 0.8),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "left") +
  scale_fill_manual(values = c("#6388b4", "#ef6f6a"), guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0,0))) #+
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 24, yend = 24, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 28, yend = 28, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 35, yend = 35, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2) +
  #annotate("segment", x = 0, xend = 40, y = 15, yend = 15, linewidth = 0.2)
  


# node search
find_sts <- as_tibble(Candida_raxml_midpoint)
find_sts <- find_sts %>% left_join(Candida_metadata %>% select(label,ST))
st_clades <- read_xlsx(metadata_file, sheet = 3)

Candida_tree <- Candida_raxml_midpoint %>% full_join(Candida_metadata, by = "label")

res_tree <- ggtree(Candida_tree,aes(size = (as.numeric(label) <95 | is.na(as.numeric(label)))))  +
  scale_size_manual(values = c(1, 0.2), guide = "none") +
  #geom_tippoint(size = 1.2, aes(color = res_cat)) + 
  geom_tiplab(size = 2, aes(label = mec_isolate_code, color = res_cat), align = TRUE, linetype = "dotted", linesize = 0.2) + 
  geom_rootedge() +
  scale_color_manual(values = c("#32a251", "#b85a0d", "#ffd94a", "grey30") , 
                     labels = c("FLC-R", "MCF-R", "MDR", "S")
                     ) +
  geom_treescale(x=0, y=0) +
  theme(legend.position = "inside",
        legend.position.inside =c(0.09,0.87), 
        legend.title = element_blank(),
        legend.margin = margin(t =0, b = 0, unit = "cm"),
        legend.text = element_text(margin = margin(l = 0, unit = "cm"))) +
  geom_cladelab(data = st_clades, 
                mapping = aes(node = st_nodes, label = st_ID),
                align = TRUE,
                offset = 0.015,
                fontsize = 2.5
  ) 

ggsave(paste0(Sys.Date(),"_",species,"_midpoint_resistance.pdf"),
       res_tree, bg="white", width = 6, height = 7.5, units = "in")

res_tree + cnvs

