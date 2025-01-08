### ---------------------------
## Purpose: plot all relative copy number data for a group of samples on a single linear ideogram
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
# Load packages
library(readxl)
library(tidyverse)
library(paletteer)
library(patchwork)

spreadsheet_list <- "~/umn/data/metadata/Calbicans_snp_depth_paths.txt"
feature_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Calbicans_SC5314_A21_plotting_features.txt"
label_file <-  "~/umn/Candida_genome_visualization/ref_genome_files/Calbicans_SC5314_A21_chr_labels.txt"
patient_data <- "~/umn/data/metadata/2024_Calbicans_sorted_patient_info.xlsx"

patient_series <- list(c("MEC279", "MEC293"), 
                       c("MEC318", "MEC319", "MEC321", "MEC322","MEC320") ,
                        c("MEC131", "MEC132", "MEC133", "MEC134", "MEC324"),
                       c("MEC352"),
                       c("MEC218", "MEC219"),
                       c("MEC157"),
                       c("MEC324"),
                       c("MEC185"),
                       c("MEC135")
                       )

patient_series <- list(c("MEC198", "MEC199", "MEC200", "MEC201"))
#patient_series <- list(c("MEC246", "MEC248", "MEC249", "MEC254"))

save_dir <-"~/umn/images/Calbicans/"

ploidy_multiplier <- 2  # this multiplied by ploidy sets the max-y scale
y_axis_labels <- c(1,2,3,4)  # manual y-axis labels, adjust as needed

inter_chr_spacing <- 150000 # size of space between chrs
chrom_outline_color <- "gray15"  # color of chromosome outlines
chrom_line_width <- 0.5  # line width of chromosome outlines

################################################################################
# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())

# Read in patient metadata
pop_data <- read_xlsx(patient_data)
pop_data[pop_data=="NA"] <- NA

# Read in and combine depth data for all isolates
depth_files <- scan(spreadsheet_list, what=character())

genome_depth <- read_xlsx(depth_files[1]) %>%
  select(chr, index, pos, plot_pos, relative_depth)
names(genome_depth)[names(genome_depth)=="relative_depth"] <-str_extract(depth_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:100){
  new_depth <- read_xlsx(depth_files[i]) %>%
    select(chr, index, pos, plot_pos, relative_depth)
  names(new_depth)[names(new_depth)=="relative_depth"] <-str_extract(depth_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_depth <- genome_depth %>%
    left_join(new_depth, by = join_by(chr,index,pos, plot_pos))
}


genome_depth <- genome_depth %>%
  pivot_longer(names_to = "sample", values_to = "relative_depth", cols=-c(chr,index,pos,plot_pos)) %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2"))

genome_depth <- genome_depth%>%
  mutate(ploidy = case_when(sample=="MEC324" ~ 3,
                                 sample=="MEC185" ~ 4,
                                 .default = 2))

genome_depth <- genome_depth %>% 
  mutate(copy_number = relative_depth * ploidy)

genome_depth <- genome_depth %>% inner_join(pop_data, by=join_by(sample))

# Small dataframes for chrom. outlines and features
chroms <- genome_depth %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  left_join(chroms, by=join_by(index))

features <- features %>%
  mutate(plot_start = start + xmin, plot_end = end + xmin)

# Tick marks to center chromosome ID label
ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs =
                  0.5, na.remove = TRUE)

################################################################################
# Plot subset of isolates and combine with patchwork 

for(k in 1:length(patient_series)){
  for(j in 1:length(patient_series[[k]])){
      tmp <- genome_depth %>%
        filter(sample==patient_series[[k]][j]) %>% 
        ggplot() +
        geom_segment(aes(x = plot_pos,
                        y = ifelse(copy_number <= 2*ploidy_multiplier, copy_number, Inf),
                         xend = plot_pos, 
                        yend = 2), 
                    color = "grey40",
                    alpha = 0.9, 
                    show.legend = FALSE) +
        geom_point(data = features, size = 1.3,
                  aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
                  position = position_nudge(y=0.1)) +
        scale_fill_manual(values = c("white", "grey26", "deepskyblue"), guide = "none") + 
        scale_shape_manual(values = c(24,21,22), guide = "none") +
        ylab(pop_data$mec_isolate_code[pop_data$sample==patient_series[[k]][j]]) +
        geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  linewidth = chrom_line_width, fill = NA,
                  colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
        scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids, position = "top") +
        scale_y_continuous(limits = c(0, 2*ploidy_multiplier), breaks = y_axis_labels) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(color = "grey90"),
              axis.ticks = element_line(color = NA),
              axis.line = element_blank(),
              axis.title.y = element_text(size = 11, color = "black", 
                                          margin = margin(r = 5, l=10),
                                          angle = 0,
                                          vjust= 0.6),
              axis.text.y = element_text(size = 8, color = "black"),
              axis.text.x = if(j==1){element_text(size=8, color ="black")
              } else{element_blank()})
  
      assign(patient_series[[k]][j], tmp)
    }
}

# Stack plots

# Triploid and tetraploid MECs
ploidy_changes <- (MEC324 + theme(plot.margin = margin(b=8)))/
  (MEC185 + theme(plot.margin = margin(b=0)))

ggsave(paste0(save_dir, Sys.Date(), "_Calbicans_MEC_ploidy_changes.pdf"),
       ploidy_changes,
       width = 6,
       height = 2,
       units = "in")

# Serial whole chromosome copy number change
aneuploids <- 
  (MEC279+ theme(plot.margin = margin(b=0)))/
  (MEC293 + theme(plot.margin = margin(b=8)))/
  (MEC318+ theme(plot.margin = margin(b=0)))/
  (MEC319+ theme(plot.margin = margin(b=0)))/
  (MEC321+ theme(plot.margin = margin(b=0)))/
  (MEC322+ theme(plot.margin = margin(b=0)))/
  (MEC320 + theme(plot.margin = margin(b=0)))
  
ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_serial_aneuploids.pdf"),
       aneuploids,
       width = 6,
       height = 5.1,
       units = "in")

# Segmental chromosome amplification
seg_aneuploid <- (MEC218+ theme(plot.margin = margin(b=0)))/
  (MEC219 + theme(plot.margin = margin(b=8)))/
  (MEC352+ theme(plot.margin = margin(b=0)))

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_segmental_aneuploidy.pdf"),
       seg_aneuploid,
       width = 6,
       height = 2.7,
       units = "in")

################################################################################
# CNV subsetting for small possible NAHR in 2 series

zoom_inter_chr_spacing <- 50000
small_cnv <- genome_depth %>% 
  filter(sample %in% c("MEC198", "MEC199", "MEC200", "MEC201"), 
         index %in% c(2,6))

chr2_size <- small_cnv %>% filter(index==2) %>% summarise(max_pos = max(pos))

small_cnv <- small_cnv %>% 
  mutate(plot_pos = case_when(index==2 ~ pos, 
                              index==6 ~(pos + chr2_size$max_pos + zoom_inter_chr_spacing)))

small_chroms <- small_cnv %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

small_features <- features <- read_tsv(feature_file, show_col_types = FALSE)

small_features <- small_features %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  filter(index %in% c(2,6)) %>% 
  left_join(small_chroms, by=join_by(index)) %>%
  mutate(plot_start = start + xmin, plot_end = end + xmin)

small_ticks <- tapply(small_cnv$plot_pos, small_cnv$index, quantile, probs =
                  0.5, na.remove = TRUE)

for(k in 1:length(patient_series)){
  for(j in 1:length(patient_series[[k]])){
    tmp <- small_cnv %>%
      filter(sample==patient_series[[k]][j]) %>% 
      ggplot() +
      geom_segment(aes(x = plot_pos,
                       y = ifelse(copy_number <= 2*ploidy_multiplier, copy_number, Inf),
                       xend = plot_pos, 
                       yend = 2), 
                   color = "grey40",
                   alpha = 0.9, 
                   show.legend = FALSE) +
      geom_point(data = small_features, size = 1.3,
                 aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
                 position = position_nudge(y=0.1)) +
      scale_fill_manual(values = c("white", "grey26", "deepskyblue"), guide = "none") + 
      scale_shape_manual(values = c(24,21,22), guide = "none") +
      ylab(pop_data$mec_isolate_code[pop_data$sample==patient_series[[k]][j]]) +
      geom_rect(data=small_chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                linewidth = chrom_line_width, fill = NA,
                colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
      scale_x_continuous(name = NULL, expand = c(0, 0), breaks = small_ticks, labels=chr_ids[c(2,6)], position = "top") +
      scale_y_continuous(limits = c(0, 2*ploidy_multiplier), breaks = y_axis_labels) +
      theme_classic() +
      theme(panel.grid.major.y = element_line(color = "grey90"),
            axis.ticks = element_line(color = NA),
            axis.line = element_blank(),
            axis.title.y = element_text(size = 11, color = "black", 
                                        margin = margin(r = 5, l=10),
                                        angle = 0,
                                        vjust= 0.6),
            axis.text.y = element_text(size = 8, color = "black"),
            axis.text.x = if(j==1){element_text(size=12, color ="black")
            } else{element_blank()})
    
    assign(patient_series[[k]][j], tmp)
  }
}


small_cnv2 <- genome_depth %>% 
  filter(sample %in% c("MEC246", "MEC248", "MEC249", "MEC254"), 
         index %in% c(1,2))

chr1_size <- small_cnv2 %>% filter(index==1) %>% summarise(max_pos = max(pos))

small_cnv2 <- small_cnv2 %>% 
  mutate(plot_pos = case_when(index==1 ~ pos, 
                              index==2 ~(pos + chr1_size$max_pos + zoom_inter_chr_spacing)))

small_chroms2 <- small_cnv2 %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

small_features2 <- features <- read_tsv(feature_file, show_col_types = FALSE)

small_features2 <- small_features2 %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  filter(index %in% c(1,2)) %>% 
  left_join(small_chroms2, by=join_by(index)) %>%
  mutate(plot_start = start + xmin, plot_end = end + xmin)

small_ticks2 <- tapply(small_cnv2$plot_pos, small_cnv2$index, quantile, probs =
                        0.5, na.remove = TRUE)

for(k in 1:length(patient_series)){
  for(j in 1:length(patient_series[[k]])){
    tmp <- small_cnv2 %>%
      filter(sample==patient_series[[k]][j]) %>% 
      ggplot() +
      geom_segment(aes(x = plot_pos,
                       y = ifelse(copy_number <= 2*ploidy_multiplier, copy_number, Inf),
                       xend = plot_pos, 
                       yend = 2), 
                   color = "grey40",
                   alpha = 0.9, 
                   show.legend = FALSE) +
      geom_point(data = small_features2, size = 1.3,
                 aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
                 position = position_nudge(y=0.1)) +
      scale_fill_manual(values = c("white", "grey26", "deepskyblue"), guide = "none") + 
      scale_shape_manual(values = c(24,21,22), guide = "none") +
      ylab(pop_data$mec_isolate_code[pop_data$sample==patient_series[[k]][j]]) +
      geom_rect(data=small_chroms2, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                linewidth = chrom_line_width, fill = NA,
                colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
      scale_x_continuous(name = NULL, expand = c(0, 0), breaks = small_ticks2, labels=chr_ids[c(1,2)], position = "top") +
      scale_y_continuous(limits = c(0, 2*ploidy_multiplier), breaks = y_axis_labels) +
      theme_classic() +
      theme(panel.grid.major.y = element_line(color = "grey90"),
            axis.ticks = element_line(color = NA),
            axis.line = element_blank(),
            axis.title.y = element_text(size = 11, color = "black", 
                                        margin = margin(r = 5, l=10),
                                        angle = 0,
                                        vjust= 0.6),
            axis.text.y = element_text(size = 8, color = "black"),
            axis.text.x = if(j==1){element_text(size=12, color ="black")
            } else{element_blank()})
    
    assign(patient_series[[k]][j], tmp)
  }
}


small_cnvs <- (MEC218+ theme(plot.margin = margin(b=0)))/
  (MEC219+ theme(plot.margin = margin(b=12)))/
  (MEC198+ theme(plot.margin = margin(b=0)))/
  (MEC201+ theme(plot.margin = margin(b=0)))/
  (MEC200+ theme(plot.margin = margin(b=0)))/
  (MEC199 + theme(plot.margin = margin(b=12)))/
  (MEC246+ theme(plot.margin = margin(b=0)))/
  (MEC248+ theme(plot.margin = margin(b=0)))/
  (MEC249+ theme(plot.margin = margin(b=0)))/
  (MEC254+ theme(plot.margin = margin(b=0)))

ggsave(paste0(save_dir,Sys.Date(),"_Calbicans_MEC_small_CNVs.pdf"),
       small_cnvs,
       width = 6,
       height = 6,
       units = "in")

################################################################################
# Plot linear genome for all isolates together - points colored by pt
p <- genome_depth %>%
  group_by(sample) %>%
  ggplot() +
  geom_point(aes(x = plot_pos,
                 y = copy_number,
                 color=as_factor(ploidy_change)),alpha = 0.8, size=0.3, shape=20,
             show.legend = FALSE) +
  scale_color_manual(values=patient_colors) +
  geom_point(data = features, size = 2,
             aes(group=index, x=plot_start, y=ymin, shape = Feature),# fill = Feature),
             position = position_nudge(y=0.07)) +
  #scale_fill_manual(values = c("white", "grey26", "deepskyblue")) +
  scale_shape_manual(values = c(17,19,15), guide = guide_legend(title = NULL)) +
  ylab('Relative copy number') +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids) +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12),
        legend.position = "bottom")

# Save plot
ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_all_copy_number.png"),
       p, width = 6, height = 2.4, units = "in", device = png, dpi = 300, bg = "white")
