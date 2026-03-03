###### Swine flu work ###### 
  ## Analyze H1 and H3 sequence data 

rm(list = ls())
wd <- "C:/Users/narma/OneDrive - University of North Carolina at Chapel Hill/PhD_research/Swine_flu/GISAID_work/sequences"
setwd(wd)

## Load packages 
library(tidyverse)
library(Biostrings)







## Read in aligned FASTA files

H1 <- readDNAStringSet(paste0(wd,"phylogeographic_analysis/inputs/sequence_files/HAH1_sequence_aligned.fasta",sep=""), format = "fasta")
H3 = readDNAStringSet(paste0(wd,"phylogeographic_analysis/inputs/sequence_files/HAH3_sequence_aligned.fasta",sep=""), format = "fasta")

## convert into matrices
H1_mat <- as.matrix(H1)
H3_mat <- as.matrix(H3)

## calculate gap percentages per sequence
gap_stats_H1 <- data.frame(
  Sequence = rownames(H1_mat),
  Gap_Percent = rowSums(H1_mat == "-" | H1_mat == "?") / ncol(H1_mat) * 100
)
gap_stats_H3 <- data.frame(
  Sequence = rownames(H3_mat),
  Gap_Percent = rowSums(H3_mat == "-" | H3_mat == "?") / ncol(H3_mat) * 100
)

#### Identified sequences to remove, remove these
H1_removeseq <- c("A/Swine/SD/Hutchinson/SS241028/2024(HA_H1)",
                  "A/Swine/SD/Hutchinson/SS241030/2024(HA_H1)",
                  "A/Swine/MO/Scott/SS241305/2024(HA_H1)",
                  "A/Swine/MO/Scott/SS241310/2024(HA_H1)")
H3_removeseq <- c("A/Swine/MO/Scott/SS241272/2024(HA_H3)",
                  "A/Swine/MO/Scott/SS241278/2024(HA_H3)",
                  "A/Swine/SD/McCook/SS241127/2024(HA_H3)",
                  "A/Swine/MN/Pipestone/SS243034/2024(HA_H3)")

H1_mat_clean <- H1_mat[!rownames(H1_mat) %in% H1_removeseq, ]
H3_mat_clean <- H3_mat[!rownames(H3_mat) %in% H3_removeseq, ]

## Convert matrices to DNAStringSet
H1_clean <- DNAStringSet(apply(H1_mat_clean, 1, paste, collapse = ""))
H3_clean <- DNAStringSet(apply(H3_mat_clean, 1, paste, collapse = ""))

####################################################
## Also for each- assign farm identifiers based on the project_id 
  # Get project_ID to farm_ID file
ID <- read.csv(paste0(wd,"phylogeographic_analysis/inputs/Project_to_Farm_ID_Mapping.csv",sep=""))

## Function to extract project ID with SS prefix
extract_project_id <- function(name) {
  # Extract /SSXXXXX/ including the SS prefix
  str_extract(name, "/SS\\d+/") %>% 
    str_remove_all("/")  # Removes the slashes -> "SS243331"
}

## Function to update names in a single DNAStringSet
update_stringset_names <- function(dna_stringset, ID) {
  current_names <- names(dna_stringset)
  project_ids <- extract_project_id(current_names)
  
  # Match project IDs (with SS) to farm IDs
  matches <- match(project_ids, ID$Project_ID)
  farm_ids <- ID$Farm_id[matches]
    # For FS1 and 1FS, change farm_id to F1S
  farm_ids <- ifelse(farm_ids %in% c("FS1", "1FS"), "F1S", farm_ids)
  
  # Create new names by appending farm_id to the end if available
  new_names <- ifelse(
    !is.na(farm_ids),
    paste0(current_names, "/", farm_ids),  # Append farm ID
    current_names  # Keep original if no match
  )
  
  names(dna_stringset) <- new_names
  return(dna_stringset)
}

## Main processor for multiple DNAStringSets
process_multiple_stringsets <- function(stringset_list, ID, output_names = NULL) {
  # Process each DNAStringSet
  updated_sets <- lapply(stringset_list, function(dss) {
    update_stringset_names(dss, ID)
  })
  
  # Apply names if provided
  if (!is.null(output_names)) {
    names(updated_sets) <- output_names
  }
  
  return(updated_sets)
}

# Prepare the list of DNAStringSets
dna_sets <- list(H1, H3, H1_clean, H3_clean)

# Process all sets (keeping SS in project IDs)
updated_sets <- process_multiple_stringsets(
  dna_sets,
  ID,
  output_names = c("H1_final", "H3_final", "H1_c_final", "H3_c_final")
)

# Extract to global environment
list2env(updated_sets, envir = .GlobalEnv)

## Export files out
writeXStringSet(H1_final, paste0(wd,"phylogeographic_analysis/outputs/HAH1_sequence_aligned_final.fasta",sep=""))
writeXStringSet(H3_final, paste0(wd,"phylogeographic_analysis/outputs/HAH3_sequence_aligned_final.fasta",sep=""))
writeXStringSet(H1_c_final, paste0(wd,"phylogeographic_analysis/outputs/HAH1_sequence_alignedsub_final.fasta",sep=""))
writeXStringSet(H3_c_final, paste0(wd,"phylogeographic_analysis/outputs/HAH3_sequence_alignedsub_final.fasta",sep=""))
