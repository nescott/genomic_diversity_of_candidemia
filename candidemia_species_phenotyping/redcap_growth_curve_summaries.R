## ---------------------------
## Purpose: Plot growth curve summaries and test correlations for all species
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

source("redcap_MIC_summary.R")

library(patchwork)
library(ggbeeswarm)
library(correlation)

save_path <- "~/umn/images/2023_growth_curves/summaries/"

growth_curves <- '58045'

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", ">32", "64", "128", "160", "256")
mcf_amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1", "2")

# function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
  formData <- list("token"=Sys.getenv("redcap_api_key"),
                   content='report',
                   format='csv',
                   report_id=report_number,
                   csvDelimiter='',
                   rawOrLabel='label',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='true',
                   returnFormat='csv'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}
################################################################################
# Read in GC report and join to MIC data

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA", 
           !primary_id %in% c("AMS5122", "AMS5123")) %>%
    select(primary_id, redcap_repeat_instance, 
           gc_date, drug_used, 
           gc_temp, gc_time, 
           k, r, t_gen, auc_l, auc_e) 

gc <- gc %>% 
    group_by(primary_id) %>% 
    filter(drug_used=="None", gc_date==max(gc_date), gc_time < 24.1)  

gc <- gc %>% 
  left_join(mic_info, by=join_by(primary_id))

# Quick summary statistics by species
gc_summary <- gc %>% 
  select(primary_id, genus_species, k, r, auc_e) %>% 
  split(.$genus_species) %>% 
  map(summary)

# For ordering species and colors for consistency
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

gc$genus_species <- factor(gc$genus_species, levels = species_count$genus_species)

################################################################################

# Add count of samples to box plots at bottom of plot
give.n <- function(x){
  return(c(y = 0, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# YPAD carrying capacity
ypad_k_metric <- ggplot(gc, aes(x=genus_species, y=k, color = genus_species)) + 
    #geom_beeswarm(cex = 0.7) +
    geom_boxplot() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 2.0)) +
    scale_color_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab(NULL) +
    ylab("OD600")

# Save YPAD carrying capacity plot
ggsave(paste0(save_path,Sys.Date(),"_MEC_YPAD_carrying_capacity.pdf"),
       ypad_k_metric,
       bg="white", 
       width = 6, 
       height= 4, 
       units="in")

# YPAD growth rate
r_metric <- ggplot(gc, aes(x=genus_species, y=r, color = genus_species)) + 
  #geom_beeswarm(size = 0.8, cex = 0.5) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values=species_colors, guide = "none") +
  theme(axis.title = element_text(color = "black",
                                  size = 12)) +
  theme(axis.text.x = element_text(angle = 30, 
                                    face="italic", 
                                    hjust = 1, 
                                    vjust = 1, 
                                    color = "black",
                                    size = 11)) +
  xlab(NULL) +
  ylab(expression(paste("Growth rate, ","hr"^-1)))

r_metric <- r_metric + stat_summary(fun.data = give.n, geom = "text", fun.y = median)

# Save growth rate plot
ggsave(paste0(save_path,Sys.Date(),"_MEC_YPAD_growth_rate.pdf"),
       r_metric,
       bg="white", 
       width = 6, 
       height= 4.2, 
       units="in")

################################################################################
# Plot both carrying capacities
k_scatter_plot <- gc %>% 
  filter(drug=="fluconazole") %>% 
  ggplot(aes(x=mean_no_drug_stationary_k, y=k, color = genus_species)) +
  geom_point(size = 0.5) +
  facet_wrap(~genus_species, nrow = 3) +
  scale_color_manual(values=species_colors, guide = "none") +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", size =9)) +
  theme(axis.text = element_text(size = 7, colour = "black")) +
  xlab("Carrying capacity in RPMI") +
  ylab("Carrying capacity in YPAD") 


ggsave(paste0(save_path,Sys.Date(),"_MEC_k_scatter.pdf"), 
       k_scatter_plot, 
       width = 6, 
       height = 5, 
       units = "in")

################################################################################
# Test for correlations
gc_corrs <- gc %>% 
  group_by(drug, genus_species) %>% 
  select(primary_id, genus_species, drug, k, mean_smg, mean_no_drug_stationary_k) %>% 
  correlation(method = "pearson")

k_vs_stationary_corrs <- gc %>% 
  group_by(drug, genus_species) %>% 
  select(primary_id, genus_species, drug, k, mean_no_drug_stationary_k) %>% 
  correlation(method = "pearson")

# Get only relevant MIC levels per drug
flc_gc <- gc %>% 
  filter(drug=="fluconazole", !primary_id %in% carrying_cap$primary_id)
flc_gc$mic50 <- factor(flc_gc$mic50, levels = flc_levels)
flc_gc$eucast_breakpoint <- factor(flc_gc$eucast_breakpoint, levels = unique(flc_gc$eucast_breakpoint))

mcf_gc <- gc %>% 
  filter(drug=="micafungin", !primary_id %in% carrying_cap$primary_id)
mcf_gc$mic50 <- factor(mcf_gc$mic50, levels = mcf_amb_levels)
mcf_gc$eucast_breakpoint <- factor(mcf_gc$eucast_breakpoint, levels = unique(mcf_gc$eucast_breakpoint))

amb_gc <- gc %>% 
  filter(drug == "amphotericin B", !primary_id %in% carrying_cap$primary_id)
amb_gc$mic50 <- factor(amb_gc$mic50, levels = mcf_amb_levels)

# Correlation for each drug
flc_gc_vs_mic <- flc_gc %>% 
  group_by(genus_species) %>% 
  mutate(mic_int = as.integer(mic50), eucast_int = as.integer(eucast_breakpoint)) %>% 
  select(genus_species, primary_id, r, mic_int, mean_no_drug_stationary_k,eucast_int ) %>% 
  correlation(method = "spearman")

mcf_gc_vs_mic <- mcf_gc %>% 
  group_by(genus_species) %>% 
  mutate(mic_int = as.integer(mic50), eucast_int = as.integer(eucast_breakpoint)) %>% 
  select(genus_species, primary_id, r, mic_int, mean_no_drug_stationary_k, eucast_int) %>% 
  correlation(method = "spearman")

amb_gc_vs_mic <- amb_gc %>% 
  group_by(genus_species) %>% 
  mutate(mic_int = as.integer(mic50)) %>% 
  select(genus_species, primary_id, r, mic_int, mean_no_drug_stationary_k) %>% 
  correlation(method = "spearman")