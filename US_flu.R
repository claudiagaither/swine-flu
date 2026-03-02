## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!

#libraries
library(ape)
library(Biostrings)
library(dplyr)
library(ggtree)
library(lubridate)
library(msa)
library(Polychrome)
library(ShortRead)
library(stringr)
library(tidyverse)
library(treeio)

## align sequences downloaded from GISAID
# Read your FASTA file
#seqs <- readDNAStringSet("C:/Users/cgait/OneDrive/Desktop/swine flu/US flu/gisaid_epiflu_sequence.fasta")
# Check for duplicate names
#dup_names <- names(seqs)[duplicated(names(seqs))]
#cat("Found duplicates:", dup_names, "\n")
# Make names unique by appending suffixes to duplicates
#names(seqs) <- make.unique(names(seqs))
#dup_names <- names(seqs)[duplicated(names(seqs))]
# Align with MAFFT (or use method = "ClustalW"/"Muscle")
#alignment <- msa(seqs, method = "ClustalW")
# Export to FASTA
#aligned_seqs <- as(alignment, "DNAStringSet")
#writeXStringSet(aligned_seqs, "aligned_CM.fasta")

## import aligned fasta to BEAUTI and export .xml file for BEAST
## alternatively do the alignment on longleaf because it's 14000 taxa !

## pull sequence dates from the aligned fasta, as BEAUTi is struggling to parse these dates?
fasta_file <- "C:/Users/cgait/OneDrive/Desktop/swine flu/US flu data/aligned_mafft.fasta"
# Read lines
lines <- readLines(fasta_file)
# Extract header lines only
headers <- lines[startsWith(lines, ">")]
# Strip the leading ">"
seq_names <- sub("^>", "", headers)
# Extract dates using regex (YYYY-MM-DD pattern) 
dates <- str_extract(seq_names, "\\d{4}-\\d{2}-\\d{2}")
# Build summary table 
date_df <- data.frame(sequence_name = seq_names, date = dates, stringsAsFactors = FALSE)
# Check for any sequences where date couldn't be parsed
#failed <- date_df[is.na(date_df$date), ]
# export to tsv to import to BEAUTi
write.table(date_df, file = "alignted_mafft_dates.tsv", sep = "\t",
  row.names = FALSE, col.names = TRUE, quote = FALSE)  


