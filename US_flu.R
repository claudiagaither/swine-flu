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

## part one: BEAST input ----

## remove problematic sequence: A/swine/Illinois/A00857131/2011|EPI_ISL_121898|A_/_H3N2||||2011-09-24|HA|4 from fasta before alignment
# Read the FASTA
#seqs <- readDNAStringSet("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v7/gisaid_epiflu_sequence.fasta")
# Check it's in there
#grep("A00857131", names(seqs), value = TRUE)
# Remove it
#seqs_clean <- seqs[!grepl("A00857131", names(seqs))]
# Verify
#length(seqs) - length(seqs_clean)  # should be 1
#length(seqs_clean)
# Write out
#writeXStringSet(seqs_clean, 
#                "C:/Users/cgait/OneDrive/Desktop/swine flu/US flu/US_flu_2010/epiflu_HA_2010_clean.fasta")


## align sequences downloaded from GISAID on longleaf using align.flu.sh

## pull and reattach sequence dates from the aligned fasta, as BEAUTi is struggling to parse
#fasta_file <- "C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v7/aligned_epiflu_HA2010.fasta"
# Read lines
#lines <- readLines(fasta_file)
# Extract header lines only
#headers <- lines[startsWith(lines, ">")]
# Strip the leading ">"
#seq_names <- sub("^>", "", headers)
# Extract dates using regex (YYYY-MM-DD pattern) 
#dates <- str_extract(seq_names, "\\d{4}-\\d{2}-\\d{2}")
# Build summary table 
#date_df <- data.frame(sequence_name = seq_names, date = dates, stringsAsFactors = FALSE)
# Check for any sequences where date couldn't be parsed
#failed <- date_df[is.na(date_df$date), ]
# export to tsv to import to BEAUTi
#write.table(date_df, file = "aligned_epiflu_HA2010_dates.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  


## N-mask non-ACTGN characters in a FASTA alignment 
#output_file <- "C:/Users/cgait/OneDrive/Desktop/nmask_aligned_epiflu_HA_2010.fasta"
#lines <- readLines(fasta_file)
#masked_lines <- character(length(lines))

#for (i in seq_along(lines)) {
#  if (startsWith(lines[i], ">")) {
#    # Header line — keep as-is
#    masked_lines[i] <- lines[i]
#  } else {
#    # Sequence line — uppercase first, then replace anything not ACTGN with N
#    seq <- toupper(lines[i])
#    seq <- gsub("[^ACTGN]", "N", seq)
#    masked_lines[i] <- seq
#  }
#}

#writeLines(masked_lines, output_file)

# Quick summary of what was changed
#original_seqs <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
#masked_seqs   <- paste(masked_lines[!startsWith(masked_lines, ">")], collapse = "")

#n_original <- nchar(gsub("[ACTGN]", "", original_seqs))
#cat("Total characters replaced with N:", n_original, "\n")



## part two: BEAST output, MCC tree ----

## taxa/sequence/tip metadata
tip_dates <- read.table("C:/Users/cgait/OneDrive/Desktop/swine flu/US flu/US_flu_2010/aligned_epiflu_HA_2010_dates.txt",
                        header = FALSE, skip = 1, col.names = c("label", "date"))
#tip_dates <- read.table("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v7/aligned_epiflu_HA2010_dates.txt",
#                        header = FALSE, skip = 1, col.names = c("label", "date"))
tip_dates$date <- as.Date(tip_dates$date)
tip_meta <- tip_dates

# Extract state name from label column (text between 2nd and 3rd "/")
tip_meta$state <- sapply(strsplit(tip_meta$label, "/"), function(x) x[3])

## clean state names
tip_meta <- tip_meta %>% mutate(state = case_when(
                                state == "A01377457" ~ "Kansas", 
                                state == "Dakota_Superior" ~ "North Dakota",
                                state == "Ilinois" ~ "Illinois",
                                state == "NE" ~ "Nebraska",
                                state == "North_Caroiina" ~ "North Carolina",
                                state == "North_Carolina" ~ "North Carolina",
                                state == "NorthCarolina" ~ "North Carolina",
                                state == "North_Dakota" ~ "North Dakota",
                                state == "NY" ~ "New York",
                                state == "South_Carolina" ~ "South Carolina",
                                state == "South_Dakota" ~ "South Dakota",
                                state == "SouthDakota" ~ "South Dakota",
                                state == "United_State" ~ NA,
                                state == "USA" ~ NA,
                                TRUE ~ state))
#table(tip_meta$state)

#broader region variable
tip_meta <- tip_meta %>% mutate(region = case_when(
    state %in% c("Alabama", "Florida", "North Carolina",
                 "South Carolina", "Tennessee", "Virginia") ~ "Southeast",
    
    state %in% c("Arizona","California",
                 "Colorado","Oklahoma", "Texas") ~ "Southwest",
    
    state %in% c("Arkansas", "Illinois", "Indiana", "Iowa","Kansas", "Kentucky", 
                 "Missouri", "Ohio","Nebraska") ~ "Midwest",
    
    state %in% c("Minnesota", "Montana","North Dakota", "South Dakota", "Wisconsin",
                 "Wyoming", "Michigan") ~ "Central North",
    
    state %in% c("Maryland", "New York", "Pennsylvania") ~ "Northeast",
    is.na(state) ~ NA_character_))
#table(tip_meta$region, useNA="always")


## once ESS values (excluding coalescent, join & prior if necessary, but ideally all values) are above 200, 
## export and combine .trees files using LogCombiner, then put combined tree into TreeAnnotator for the MCC tree
## resampling every 100,000 states instead of 10,000 in the original run to thin 90%, 10,000 was very slow


## import combined maximum clade credibility (MCC) tree from TreeAnnotator
#mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v7/mcc_comb_epiflu_HA2010_v7_10k.trees")
mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v6/mcc_comb_epiflu_HA2010_v6_50k.trees")
options(ignore.negative.edge = TRUE)

## import combined log file for clock rates
log <- read.table("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v6/comb_epiflu_HA2010_v6.log",
                  header = TRUE, comment.char = "#", sep = "\t")

#names(log)  # look for clock.rate or clockRate
#mean(log$clock.rate)  # or median — should be ~0.002–0.005 for flu HA

# identical rates for strict clock, different for relaxed clock
#td <- as_tibble(mcc_tree)
#summary(td$rate)
# should have specified a strict clock but we can double check? 
# or maybe they are pretty much the same
# and the variation is due to exponential growth coalescent? but they are all super close to 0.004

# extract and view data
tree_phylo <- mcc_tree@phylo

## Posterior labels (only show > 0.999)
#mcc_tree@data$posterior_label <- ifelse(mcc_tree@data$posterior > 0.999, round(mcc_tree@data$posterior, 3), NA)

# Each row in td$rate corresponds to that node's parent branch
# The edge table maps [parent, child] — index by child node
#tree_phylo$edge.length

# Find tips on the longest terminal branches
#td <- as_tibble(mcc_tree)
#tip_rows <- td %>% filter(!is.na(label)) %>% 
#  arrange(desc(branch.length)) %>% 
#  head(20)
#print(tip_rows %>% select(label, branch.length))

# Check — should be ~133
#max(node.depth.edgelength(mcc_tree@phylo))

# color by states for initial sanity check
tree_states <- ggtree(mcc_tree) %<+% tip_meta +
  geom_tippoint(aes(color = state), size = 2) +
  theme_tree2() + coord_cartesian(xlim = c(105, 135)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 2))
#tree_states

#color tips by broader region
tree_region <- ggtree(mcc_tree) %<+% tip_meta +
  geom_tippoint(aes(color = region), alpha=0.5, size = 3) +
  scale_color_brewer(palette="Accent") +
  theme_tree2() + coord_cartesian(xlim = c(0, 135)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) +
  guides(color = guide_legend(ncol = 1))
tree_region
