## ---------------------------
## Purpose: Summarize and plot MEC isolate (EUCAST) MIC and SMG data by species 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------

source("redcap_MIC_summary.R")

## Load packages
library(patchwork)
library(ggbeeswarm)

save_path <- "~/umn/images/2023_MICs/summaries/"

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))

mic_info <- mic_info %>% 
  mutate(mic50 = case_when((drug=="micafungin" & genus_species=="C. parapsilosis" & mic50=="2") ~ ">1", 
                           .default = mic50))


flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", ">32")
mcf_amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1")

mic_flc <- mic_info %>% 
  filter(drug=="fluconazole")
mic_flc$mic50 <- factor(mic_flc$mic50, levels = flc_levels)

mic_mcf <- mic_info %>% 
  filter(drug=="micafungin")
mic_mcf$mic50 <- factor(mic_mcf$mic50, levels=mcf_amb_levels)

mic_amb <- mic_info %>% 
  filter(drug=="amphotericin B")
mic_amb$mic90 <- factor(mic_amb$mic90, levels = mcf_amb_levels)

################################################################################
# Make 3 MIC plots - subsets of species as columns, drugs as rows (patchwork)
flc1 <- ggplot(mic_flc %>% filter(genus_species %in% c("C. albicans", 
                                                       "C. glabrata", 
                                                       "C. parapsilosis",
                                                       "C. lusitaniae") ), aes(x=mic50)) + 
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(NULL)+
  ylab("\n\nFluconazole") +
  theme(strip.text = element_text(face = "italic", size =10))+
  theme(panel.grid.major.x = element_blank()) + 
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4, size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  geom_vline(data=filter(mic_flc, genus_species %in% c("C. albicans", 
                                                        "C. lusitaniae", 
                                                        "C. parapsilosis"
                                                        )), 
              aes(xintercept = "4"), linetype =2) +
  geom_vline(data = filter(mic_info, genus_species=="C. glabrata"),
              aes(xintercept = "16"), linetype = 2)

mcf1 <- ggplot(mic_mcf %>% filter(genus_species %in% c("C. albicans", 
                                                       "C. glabrata", 
                                                       "C. parapsilosis",
                                                       "C. lusitaniae")), aes(x=mic50)) + 
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  facet_wrap(.~genus_species, nrow = 1, scales = "free_x")+
  theme_bw() +
  xlab(NULL)+
  ylab("Number of isolates\n\nMicafungin")+
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) + 
  geom_vline(data=filter(mic_mcf, genus_species=="C. albicans"), 
             aes(xintercept = "0.016"),linetype = 2) +
  geom_vline(data = filter(mic_mcf, genus_species=="C. glabrata"),
             aes(xintercept = "0.032"), linetype = 2) +
  geom_vline(data = filter(mic_mcf, genus_species=="C. parapsilosis"),
             aes(xintercept = ">1"), linetype = 2) 

amb1 <- ggplot(mic_amb %>% filter(genus_species %in% c("C. albicans", 
                                                       "C. glabrata", 
                                                       "C. parapsilosis",
                                                       "C. lusitaniae")), aes(x=mic90)) + 
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(expression(paste("MIC, ", mu,"g/mL"))) +
  ylab("\n\nAmphotericin B") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle=90, vjust=-0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) + 
  geom_vline(data=filter(mic_amb, genus_species %in% c("C. albicans", 
                                                        "C. glabrata", 
                                                        "C. parapsilosis" 
                                                        )), 
              aes(xintercept = "1"), linetype =2)

full_mic_plot1 <- flc1/mcf1/amb1

ggsave(paste0(save_path,Sys.Date(),"_MEC_MIC_summary1.pdf"),
       full_mic_plot1,
       bg="white", 
       width=5.5, 
       height = 6.5, 
       units = "in")

################################################################################
flc2 <- mic_flc %>% filter(genus_species %in% c("C. krusei", 
                                                "C. orthopsilosis", 
                                                "C. kefyr",
                                                "C. tropicalis")) %>% 
                           ggplot(aes(x=mic50)) +
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = flc_levels) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  facet_wrap(.~genus_species, nrow = 1) +
  theme_bw() +
  xlab(NULL)+
  ylab("\n\nFluconazole") +
  theme(strip.text = element_text(face = "italic", size =10))+
  theme(panel.grid.major.x = element_blank()) + 
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4, size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  geom_vline(data=filter(mic_flc, genus_species %in% c("C. tropicalis", 
                                                        "C. kefyr", 
                                                        "C. orthopsilosis"
                                                       )), 
             aes(xintercept = "4"), linetype =2) 


mcf2 <- ggplot(mic_mcf %>% filter(genus_species %in% c("C. krusei", 
                                                       "C. orthopsilosis", 
                                                       "C. kefyr",
                                                       "C. tropicalis")), aes(x=mic50)) + 
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(NULL)+
  ylab("Number of isolates\n\nMicafungin")+
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) 
  

amb2 <- ggplot(mic_amb %>% filter(genus_species %in% c("C. krusei", 
                                                       "C. orthopsilosis", 
                                                       "C. kefyr",
                                                       "C. tropicalis")), aes(x=mic50)) +
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(expression(paste("MIC, ", mu,"g/mL"))) +
  ylab("\n\nAmphotericin B") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle=90, vjust=-0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) + 
  geom_vline(data=filter(mic_amb, genus_species %in% c("C. krusei", 
                                                       "C. tropicalis")), 
                                                        
                                                       
             aes(xintercept = "1"), linetype =2)

full_mic_plot2 <- flc2/mcf2/amb2

ggsave(paste0(save_path,Sys.Date(),"_MEC_MIC_summary2.pdf"),
       full_mic_plot2,
       bg="white", 
       width=5.5, 
       height = 6.5, 
       units = "in")

################################################################################
flc3 <- ggplot(mic_flc %>% filter(genus_species %in% c("C. utilis", 
                                                       "C. dubliniensis", 
                                                       "C. nivariensis") ), aes(x=mic50)) +
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = flc_levels) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(NULL)+
  ylab("\n\nFluconazole") +
  theme(strip.text = element_text(face = "italic", size =10)) +
  theme(panel.grid.major.x = element_blank()) + 
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4, size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  geom_vline(data=filter(mic_flc, genus_species %in% c( "C. dubliniensis",
                                                         "C. nivariensis",
                                                         "C. utilis")), 
                                                        
             aes(xintercept = "4"), linetype =2) 

mcf3 <- ggplot(mic_mcf %>% filter(genus_species %in% c("C. utilis", 
                                                       "C. dubliniensis", 
                                                       "C. nivariensis")), aes(x=mic50)) +
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(NULL)+
  ylab("Number of isolates\n\nMicafungin")+
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = -0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank())  

amb3 <- ggplot(mic_amb %>% filter(genus_species %in% c("C. utilis", 
                                                       "C. dubliniensis", 
                                                       "C. nivariensis")), aes(x=mic90)) + 
  geom_bar(aes(fill=genus_species), just = 1) +
  scale_fill_manual(values=species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  scale_y_continuous(limits = c(0,4)) +
  facet_wrap(.~genus_species, nrow = 1)+
  theme_bw() +
  xlab(expression(paste("MIC, ", mu,"g/mL"))) +
  ylab("\n\nAmphotericin B") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size=11, color = "black")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle=90, vjust=-0.4,size = 9, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) + 
  geom_vline(data=filter(mic_amb, 
                         genus_species=="C. dubliniensis"), 
             aes(xintercept = "1"), linetype =2)

full_mic_plot3 <- flc3/mcf3/amb3

ggsave(paste0(save_path,Sys.Date(),"_MEC_MIC_summary3.pdf"),
       full_mic_plot3,
       bg="white", 
       width=4.5, 
       height = 6.5, 
       units = "in")
################################################################################
# Add count of samples to box plots at bottom of plot
give.n <- function(x){
  return(c(y = 0, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


smg_plot <- mic_info %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug =="fluconazole") %>% 
  ggplot(aes(x=genus_species, y=mean_smg, color = genus_species)) + 
  geom_boxplot() +  
  #geom_beeswarm(size=1, cex=0.8) +#geom_violin() +
  #geom_point(data = filter(mic_info, genus_species=="C. nivariensis")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_color_manual(values=species_colors, guide = "none") +
  theme(axis.text.x = element_text(angle = 30, 
                                    face="italic", 
                                    hjust = 1, 
                                    vjust = 1,
                                   color = "black",
                                   size = 10)) +
  xlab(NULL) +
  ylab("Fluconazole SMG")

smg_plot <-smg_plot + stat_summary(fun.data = give.n, geom = "text", fun.y = median)

ggsave(paste0(save_path, Sys.Date(),"_MEC_FLC_SMG_summary.pdf"),
       smg_plot,
       bg="white", 
       width = 6, 
       height=4.2, 
       units="in")

###############################################################################
no_drug_k <- ggplot(control_od, aes(x=genus_species, y=mean_stationary_k, color = genus_species)) +
  geom_boxplot() +
  theme_bw() +
  #scale_y_continuous(limits = c(0)) +
  scale_color_manual(values=species_colors, guide = "none") +
  theme(axis.text.x = element_text(angle = 30, 
                                   face="italic", 
                                   hjust = 1, 
                                   vjust = 1,
                                   colour = "black",
                                   size = 10)) +
  xlab(NULL) +
  ylab("OD530")

no_drug_k <- no_drug_k + stat_summary(fun.data = give.n, geom = "text", fun.y = median)
ggsave(paste0(save_path, Sys.Date(),"_MEC_no_drug_carrying_capacity.pdf"),
       no_drug_k,
       bg="white", 
       width = 6 , 
       height=4.2, 
       units="in")
