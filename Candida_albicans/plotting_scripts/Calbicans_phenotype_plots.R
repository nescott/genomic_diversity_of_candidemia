## ---------------------------
## Purpose: Summary plots of growth curve metrics
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

source("ch3/Calbicans_redcap_summary.R")
library(patchwork)
library(ggbeeswarm)

calbicans_isolates <- sample_info %>% 
  filter(genus_species=="C. albicans") %>%
  select(primary_id, genus_species)

drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))

# Growth rate
gc_r <- gc %>% right_join(calbicans_isolates) %>% 
  ggplot(aes(x = genus_species, y = r)) +
  geom_beeswarm(size=1, cex = 3, color = "#88CCEE") +
  ylab(expression(paste("Growth rate, ","hr"^-1))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Carrying capacity in YPAD
gc_k <- gc %>% right_join(calbicans_isolates) %>% 
  ggplot(aes(x = genus_species, y = k)) +
  geom_beeswarm(size=1, cex = 3, color = "#72b7f4") +
  ylab("OD600, YPAD") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Carrying capacity in RPMI
rpmi_k <- control_od %>% filter(genus_species=="C. albicans") %>% 
  ggplot(aes(x = genus_species, y = mean_stationary_k)) +
  geom_beeswarm(size=1, cex = 3, color = "#5ba1f9") +
  ylab("OD530, RPMI") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

flc_smg <- ggplot(mic_info %>% filter(genus_species == "C. albicans", drug=="fluconazole"), 
                  aes(x = genus_species, y = mean_smg)) +
 geom_beeswarm(size=1, cex = 3, color = "#458CFF") +
  ylab("Mean supra-MIC growth") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

all_phenos <- gc_r + gc_k + rpmi_k + flc_smg + plot_annotation(tag_levels = "A")

ggsave(paste0(Sys.Date(),"_Calbicans_MEC_phenotypes.pdf"),
       all_phenos,
       width = 6,
       height = 4,
       units = "in")

###############################################################################
# MIC sub-plots
flc <- ggplot(mic_info %>% filter(genus_species == "C. albicans", drug=="fluconazole"), aes(x=mic50))+
    geom_bar(fill = "#88CCEE", just = 1) +
    #scale_fill_manual(values = species_colors, guide = "none") +
    scale_x_discrete(limits = c("0.5", "1", "2", "4", "8", "16", "32", ">32")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1.4, size = 8)) +
    xlab("\nFluconazole MIC50") +
    ylab(NULL)

mcf <- ggplot(mic_info %>% filter(genus_species == "C. albicans", drug=="micafungin"), aes(x=mic50))+
    geom_bar(fill = "#67acf7", just = 1) +
    #scale_fill_manual(values = species_colors, guide = "none") +
    scale_x_discrete(limits = c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1"),
                     labels = c ("0.01", "0.03", "0.06", "0.125", "0.25 ", "0.5 ", "1  ", ">1 ")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, size = 8)) +
    xlab("\nMicafungin MIC50") +
    ylab(NULL)

amb <- ggplot(mic_info %>% filter(genus_species == "C. albicans", drug=="amphotericin B"), aes(x=mic50))+
    geom_bar(fill = "#458CFF", just = 1) +
    #scale_fill_manual(values = species_colors, guide = "none") +
    scale_x_discrete(limits = c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1"),
                     labels = c ("0.01", "0.03", "0.06", "0.125", "0.25 ", "0.5 ", "1  ", ">1  ")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, size = 8)) +
    xlab("\nAmphotericin B MIC90") +
    ylab(NULL)

flc + mcf + amb
ggsave("Calbicans_MEC_MICs.png", bg="white", width = 11, height = 4, units = "in", dpi=300)

#flc_smg <- ggplot(mic_info %>% filter(genus_species == "C. albicans", drug=="fluconazole"), 
#                  aes(x = genus_species, y = mean_smg)) +
#  geom_pointrange(aes(ymax=(mean_smg + sd_smg), 
#                      ymin=(mean_smg - sd_smg)), 
#                  position = position_jitter(width = 0.5, height = 0.1, seed = 123),
#                  alpha=0.6, 
#                  show.legend = FALSE) +
#  ylab("Mean supra-MIC growth") +
#  theme_bw() +
#  theme(axis.title.x = element_blank(),
#        axis.text.x = element_blank())
