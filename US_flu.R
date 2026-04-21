## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!


## part one: BEAST input ----

## part A: BEAST input
#libraries
library(ape)
library(Biostrings)
library(dplyr)
library(ggtree)
library(lubridate)
library(msa)
library(paletteer)
library(Polychrome)
library(phangorn)
library(readxl)
library(scico)
library(sf)
library(ShortRead)
library(stringr)
library(tidyverse)
library(treeio)
library(viridis)
library(writexl)


## import US climate metadata
states <- read_sf("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/US_State_Boundaries/US_State_Boundaries.shp")
## remove non continental US states
states <- states %>% filter(NAME!="District of Columbia")
states <- states %>% filter(NAME!="U.S. Virgin Islands")
states <- states %>% filter(NAME!="Puerto Rico")

## import clade assignments from BVBRC
tip_clades <- read.csv("C:/Users/cgait/OneDrive/Desktop/1990_v1/clade_assignment_updated.csv")

## import combined maximum clade credibility (MCC) tree from TreeAnnotator
mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/1990_v1/mcc_1990_v1_150k.trees")
options(ignore.negative.edge = TRUE)

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
#fasta_file <- "C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990 US flu/v1/H3_aligned_deduped_nmask_final.fasta"
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
#write.table(date_df, file = "aligned_HA_1990_dates.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  


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


## import & clean metadata
tip_dates <- read.table("C:/Users/cgait/OneDrive/Desktop/1990_v1/aligned_HA_1990_dates.txt",
                        header = TRUE, sep = "\t")
tip_dates$date <- as.Date(tip_dates$date)
tip_meta <- tip_dates

# Extract state name from label column (text between 2nd and 3rd "/")
tip_meta$state <- sapply(strsplit(tip_meta$sequence_name, "/"), function(x) x[3])

## clean state names
tip_meta <- tip_meta %>% mutate(state = case_when(
                                state == "A01377457" ~ "Kansas", 
                                state == "Dakota_Superior" ~ "North Dakota",   ## let's double check this!
                                state == "IA" ~ "Iowa",
                                state == "IN" ~ "Indiana",
                                state == "Ilinois" ~ "Illinois",
                                state == "Gainesville" ~ "Florida",
                                state == "NE" ~ "Nebraska",
                                state == "MI" ~ "Michigan",
                                state == "ISU-Missouri" ~ "Missouri",
                                state == "Missourri" ~ "Missouri",
                                state == "MN" ~ "Minnesota",
                                state == "North_Caroiina" ~ "North Carolina",
                                state == "North_Carolina" ~ "North Carolina",
                                state == "NorthCarolina" ~ "North Carolina",
                                state == "North_Dakota" ~ "North Dakota",
                                state == "New_Mexico" ~ "New Mexico",
                                state == "NY" ~ "New York",
                                state == "New_York" ~ "New York",
                                state == "OK" ~ "Oklahoma",
                                state == "Okhaloma" ~ "Oklahoma",
                                state == "South_Carolina" ~ "South Carolina",
                                state == "South_Dakota" ~ "South Dakota",
                                state == "SouthDakota" ~ "South Dakota",
                                state == "West_Virginia" ~ "West Virginia",
                                state == "Winconsin" ~ "Wisconsin",
                                state == "United_State" ~ NA,
                                state == "United State" ~ NA,
                                state == "USA" ~ NA,
                                TRUE ~ state))
#table(tip_meta$state, useNA="always")

#broader region variable
tip_meta <- tip_meta %>% mutate(region = case_when(
    state %in% c("Alabama","Florida","Georgia","North Carolina","West Virginia",
                 "South Carolina", "Tennessee", "Virginia","Louisiana") ~ "Southeast",
    
    state %in% c("Arizona","California","New Mexico","Oregon",
                 "Colorado","Oklahoma", "Texas","Utah") ~ "Southwest/West",
    
    state %in% c("Arkansas","Illinois","Indiana","Iowa","Kansas","Kentucky", 
                 "Missouri","Ohio","Nebraska") ~ "Midwest",
    
    state %in% c("Minnesota","Montana","North Dakota","South Dakota","Wisconsin",
                 "Wyoming","Michigan") ~ "Central North",
    
    state %in% c("Maryland", "New York", "Pennsylvania") ~ "Northeast",
    is.na(state) ~ NA_character_))

#table(tip_meta$region, useNA = "always")

## filter out states with under 5 sequences
#tip_meta <- tip_meta %>% mutate(state = case_when( state %in% c("Arizona","Florida","Georgia","Louisiana",
#               "New York","New Mexico","Oregon","South Carolina","Tennessee","West Virginia") ~ NA, 
#               TRUE ~ state))
#table(tip_meta$state, useNA="always")

# Create complete metadata for all tree tips
all_tree_tips <- data.frame(sequence_name = mcc_tree@phylo$tip.label)

# Merge with your existing tip_meta (keeping all tree tips)
tip_meta <- merge(all_tree_tips, tip_meta, by = "sequence_name", all.x = TRUE)

# Then merge with clade assignments
tip_meta_clades <- merge(tip_meta, tip_clades, by = "sequence_name", all.x = TRUE)

#clean clade assignments
tip_meta_clades <- tip_meta_clades %>% mutate(clade = case_when(
                                       clade == "Other-Human-2020" ~ "Other human",
                                       clade == "Other-Human-1970" ~ "Other human",
                                       clade == "Other-Human-1970-like" ~ "Other human",
                                       clade == "Other-Human-2000" ~ "Other human",
                                       clade == "Other-Human-2000-like" ~ "Other human",
                                       clade == "Other-Human-2010" ~ "Other human",
                                       clade == "Other-Human-2010-like" ~ "Other human",
                                       TRUE ~ clade))
tip_meta_clades <- tip_meta_clades %>% mutate(clade = case_when(
                                      clade == "1990.4-like" ~ "1990.4(-like)",
                                      clade == "1990.4" ~ "1990.4(-like)",
                                      clade == "1990.4.b1" ~ "1990.4.b1 & b2",
                                      clade == "1990.4.b2" ~ "1990.4.b1 & b2",
                                      clade == "" ~ "Missing",
                                      clade == "2010.1-like" ~ "2010.1(-like)",
                                      clade == "2010.1" ~ "2010.1(-like)",
                                      is.na(clade) ~ "Missing",
                                      TRUE ~ clade))
#table(tip_meta_clades$clade, useNA="always")

## once ESS values (excluding coalescent, join & prior if necessary, but ideally all values) are above 200, 
## export and combine .trees files using LogCombiner, then put combined tree into TreeAnnotator for the MCC tree

## part two: MCC tree---- 

## prune sequences from 1970s? or just crop tree ??
#remove_tips <- c("A/swine/Illinois/A00857131/2011|EPI_ISL_121898|A_/_H3N2||||2011-09-24|HA|4")
# verify both tips actually exist in the tree before pruning
#missing <- setdiff(remove_tips, mcc_tree@phylo$tip.label)
#if (length(missing)) stop(paste("These tips are not in the tree:", paste(missing, collapse = ", ")))

## prune tree — drop.tip removes the named tips directly
#mcc_tree <- treeio::drop.tip(mcc_tree, remove_tips)
# remove dates that are not in the tree
#tip_meta <- tip_meta[tip_meta$sequence_name %in% mcc_tree@phylo$tip.label, , drop = FALSE]
# verify
#stopifnot(length(mcc_tree@phylo$tip.label) == nrow(tip_meta))
# prune metadata to match
#tip_meta <- tip_meta[!tip_meta$sequence_name %in% remove_tips, , drop = FALSE]

# extract and view data
tree_phylo <- mcc_tree@phylo

#table(tip_meta_clades$clade, useNA="always")

#color tips by clade assignments made using BV-BRC
tree_clade <- ggtree(mcc_tree) %<+% tip_meta_clades +
  geom_tippoint(aes(color = clade), alpha=0.5, size = 3) +
  scale_color_paletteer_d("ggsci::default_ucscgb") + 
  theme_tree2() + coord_cartesian(xlim = c(1990, 2033)) +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) + guides(color = guide_legend(ncol = 1))

#tree_clade

# See how many sequences match between tree and clade assignments
tree_tips <- mcc_tree@phylo$tip.label
clade_names <- tip_clades$sequence_name
length(intersect(tree_tips, clade_names))

# See which clade assignments are missing (not matching tree tips)
missing_in_tree <- setdiff(clade_names, tree_tips)

# See which tree tips don't have clade assignments
missing_clades <- setdiff(tree_tips, clade_names)

## prune missing_clades tips from tree and re-plot?
pruned_tree <- drop.tip(mcc_tree@phylo, missing_clades)

# If mcc_tree is a treedata object, update it
mcc_tree_pruned <- mcc_tree
mcc_tree_pruned@phylo <- pruned_tree

# Re-plot with pruned tree
tree_clade_pruned <- ggtree(mcc_tree_pruned) %<+% tip_meta_clades +
  geom_tippoint(aes(color = clade), alpha = 0.5, size = 3) +
  scale_color_paletteer_d("ggsci::default_ucscgb") +
  theme_tree2() + coord_cartesian(xlim = c(120, 175)) +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) +
  guides(color = guide_legend(ncol = 1))

#tree_clade_pruned
#ggsave("C:/Users/cgait/OneDrive/Desktop/clades_pruned.jpeg",width=20,height=25,units=c("cm"),tree_clade_pruned)

## part three: continuous traits? ----

## density of samples from each state
## pull table of all unique state names and number of sequences from each, including those with no state
state_counts <- as.data.frame(table(tip_meta$state, useNA = "ifany"))
colnames(state_counts) <- c("state", "n_sequences")
all_states <- data.frame(state = state.name)
state_counts_full <- merge(all_states, state_counts, by = "state", all.x = TRUE)

## merge with states shapefile
states$state <- states$NAME
states <- left_join(states, state_counts_full, by="state")
states$log_n <- log(states$n_sequences)
remove(all_states, state_counts, state_counts_full)
states <- states %>% filter(state != "Alaska")
states <- states %>% filter(state != "Hawaii")

## filter states for inclusion
states <- states %>% filter(!is.na(n_sequences))
#sum(states$n_sequences)
states <- states %>% filter(n_sequences>=10)

## hog centers for each state
## pulled just looking at the ESRI map lol will do more formally later
locations <- data.frame(state = c("Iowa", "Minnesota", "Illinois", 
                                  "Indiana", "Ohio","North Carolina", 
                                  "Pennsylvania", "Michigan", "Wisconsin",
                                  "Kentucky", "Arkansas", "Oklahoma", 
                                  "Texas", "Missouri","Kansas", 
                                  "Colorado", "Wyoming", "Utah", 
                                  "Nebraska","South Dakota"),
                        lat = c(42.78337081346157, 44.10462846785306, 40.38267144729745,
                                40.84318098374256, 40.61332295049134, 35.235921,
                                40.47488269049696, 42.95188067992797, 43.630949241417866,
                                37.23207623324495, 34.29984836061682, 36.666475839046754,
                                31.542043656011728, 40.073903437901684, 39.635912212028494,
                                39.02144609276473, 41.50905940171905, 38.17199803139623,
                                41.452249465481124, 43.37415915245506),
                        lon = c(-94.42785719132587, -95.26198230392124, -89.27144699793,
                                -86.31409432600094, -83.63478335827033, -78.192904,
                                -76.93887362132278, -85.80965089379956, -90.18386220641291,
                                -87.14821943859631, -93.87177340330743, -99.81175526572906,
                                -96.5763615952025, -92.10241780946365, -96.88224256223276,
                                -102.38996085242769, -104.46263573362431, -113.66328849073699,
                                -97.28410432335959, -97.48631647186757))
state_centers <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

US_samples <- ggplot() + geom_sf(data = states, aes(fill = n_sequences)) +
              geom_sf(data = state_centers, color="pink4", fill="gold", size=3, shape=22) +
              scale_fill_scico(palette = "hawaii", direction = -1) +
              theme_classic() 
#US_samples



## attach climate data over time ?


