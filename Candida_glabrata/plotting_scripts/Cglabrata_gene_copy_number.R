################################################################################
## Purpose: Overlap GC-corrected read depth with gene locations to find candidate
## copy number changes
## Author: Nancy Scott
## Email: scot0854@umn.edu
################################################################################

library(readxl)
library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)

# Sample inputs
spreadsheet_list <- "ch4/Cglabrata_depth_files.txt"
copy_number_files <- scan(spreadsheet_list, what=character())

# Genome inputs
candida_gff <- "~/umn/local/ref_genomes/C_glabrata_CBS138_version_s05-m03-r02_features2.gff"
region_file <- "~/umn/local/ref_genomes/C_glabrata_CBS138_version_s05-m03-r02_regions.bed"

# Function returns GRange of gene candidates, either above or below relative_depth threshold, with minimum specified overlap
coverage_function <- function(depth_range, gene_range, op, threshold, overlap_fraction){
  cov1 <-depth_range[get(op)(depth_range$relative_depth, threshold)]
  cov <- reduce(cov1)
  cov <- cov[width(cov)>500]
  
  gene_hits <- findOverlaps(genes, cov, ignore.strand = TRUE)
  gene_overlaps <- pintersect(genes[queryHits(gene_hits)], cov[subjectHits(gene_hits)])
  percent_overlap <- width(gene_overlaps) /width(genes[queryHits(gene_hits)])
  hits <- gene_hits[percent_overlap > overlap_fraction]
  
  genes_of_interest <- gene_range[queryHits(hits)]
  genes_of_interest <- genes_of_interest[genes_of_interest$orf_classification != "Dubious"]
  return(genes_of_interest)
}


# GFF
features <- import.gff(candida_gff)
genes <- features[features$type=="gene"]

# Ranges of interest (not repeat regions, LTRs, retrotransposons)
regions <- import.bed(region_file)

# Read in 
for(i in 1:length(copy_number_files)){
  sample_id <- str_extract(copy_number_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  # Make Granges of depth data
  genome_depth <- read_xlsx(copy_number_files[i]) 
  colnames(genome_depth)[2] <- "begin"

  genome_depth <- genome_depth %>% 
    mutate(end = begin + 499) %>% 
    select(chr,begin,end,relative_depth) %>% 
    mutate(relative_depth = case_when(is.na(relative_depth) ~ 0,
                                      .default = relative_depth))
  
  genome_mean_depth <- genome_depth %>% 
    group_by(chr) %>% 
    summarise(avg_depth = mean(relative_depth))
  
  genome_mean_depth_keep <- genome_mean_depth %>% 
    filter(avg_depth < 1.28)
  
  genome_cov <- toGRanges(as.data.frame(genome_depth))

  # Subset to only genes
  genome_cov <- subsetByOverlaps(genome_cov, 
                                 genes)
  overlaps <- findOverlaps(query = genes,
                           subject = genome_cov)
  
  mcols(genome_cov)$orf[subjectHits(overlaps)] <- mcols(genes)$ID[queryHits(overlaps)]
  mcols(genome_cov)$Gene[subjectHits(overlaps)] <- mcols(genes)$Gene[queryHits(overlaps)]
  
  # Subset to regions of interest
  genome_cov <- subsetByOverlaps(genome_cov, regions, 
                                minoverlap = 250)

  # Get low coverage genes of interest, make dataframe, add sample ID and write csv
  low_cov <- coverage_function(genome_cov, genes, "<", 0.25, 0.8)
  lc_df <- as.data.frame(low_cov)
  lc_df <- lc_df %>% select(-c(Note,Alias,Parent))
  lc_df$sample <- c(rep(sample_id, times = length(lc_df$seqnames)))
  write.table(lc_df, 
            "Cglabrata_MEC_candidate_deletions.tsv", 
            append = TRUE, 
            quote = FALSE,
            sep = "\t",
            col.names = case_when(i==1 ~ TRUE,
                                  .default = FALSE))
  
  avg_del_depth <- as.data.frame(subsetByOverlaps(genome_cov, 
                                                  low_cov, 
                                                  ignore.strand = TRUE, 
                                                  minoverlap = 100)) %>% 
    group_by(orf) %>% summarise(mean_rel_depth = mean(relative_depth))
  avg_del_depth$sample <- c(rep(sample_id, times=length(avg_del_depth$orf)))
  
  write.table(avg_del_depth, 
            "Cglabrata_MEC_candidate_deletion_coverage.tsv", 
            append = TRUE, 
            quote = FALSE,
            sep = "\t",
            col.names = case_when(i==1 ~ TRUE,
                                  .default = FALSE))


  # Get high coverage genes of interest, append to csv
  high_cov <- coverage_function(genome_cov, genes, ">", 1.9, 0.75)
  hc_df <- as.data.frame(high_cov)
  hc_df$sample <- c(rep(sample_id, times = length(hc_df$seqnames)))
  hc_keep <- hc_df %>% 
    filter(seqnames %in% genome_mean_depth_keep$chr) %>% 
    select(-c(Note,Alias,Parent))
  write.table(hc_keep, 
              "Cglabrata_MEC_candidate_duplications.tsv", 
             append = TRUE, 
              quote = FALSE, 
              sep = "\t",
              col.names = case_when(i==1 ~ TRUE,
                                    .default = FALSE))
  
  
  avg_dup_depth <- as.data.frame(subsetByOverlaps(genome_cov, 
                                                  high_cov, 
                                                  ignore.strand = TRUE, 
                                                  minoverlap = 100)) %>%
    filter(seqnames %in% genome_mean_depth_keep$chr) %>% 
    group_by(orf) %>% summarise(mean_rel_depth = mean(relative_depth))
  avg_dup_depth$sample <- c(rep(sample_id, times=length(avg_dup_depth$orf)))
  
  write.table(avg_dup_depth, 
              "Cglabrata_MEC_candidate_duplication_coverage.tsv", 
             append = TRUE, 
              quote = FALSE, 
              sep = "\t",
              col.names = case_when(i==1 ~ TRUE,
                                    .default = FALSE))

}
