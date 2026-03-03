## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!

#libraries
library(ape)
#library(Biostrings)
library(dplyr)
library(ggtree)
library(lubridate)
library(msa)
library(Polychrome)
library(ShortRead)
library(stringr)
library(tidyverse)
library(treeio)

## align sequences downloaded from GISAID (on longleaf as it is bulky)

## pull and reattach sequence dates from the aligned fasta, as BEAUTi is struggling to parse

## N-mask non-ACTGN characters in a FASTA alignment 

fasta_file <- "C:/Users/cgait/OneDrive/Desktop/swine flu/US flu/US flu data/aligned_mafft.fasta"       # input
output_file <- "C:/Users/cgait/OneDrive/Desktop/swine flu/US flu/US flu data/aligned_mafft_nmask.fasta"  # output

lines <- readLines(fasta_file)

masked_lines <- character(length(lines))

for (i in seq_along(lines)) {
  if (startsWith(lines[i], ">")) {
    # Header line — keep as-is
    masked_lines[i] <- lines[i]
  } else {
    # Sequence line — uppercase first, then replace anything not ACTGN with N
    seq <- toupper(lines[i])
    seq <- gsub("[^ACTGN]", "N", seq)
    masked_lines[i] <- seq
  }
}

writeLines(masked_lines, output_file)

# Quick summary of what was changed
original_seqs <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
masked_seqs   <- paste(masked_lines[!startsWith(masked_lines, ">")], collapse = "")

n_original <- nchar(gsub("[ACTGN]", "", original_seqs))
cat("Total characters replaced with N:", n_original, "\n")
