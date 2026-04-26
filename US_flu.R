## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!


## part one: baseline MCC tree ----

source("C:/Users/cgait/OneDrive/Desktop/swine flu/US_flu_functions.R")

## part one A: BEAST input
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
library(xml2)


## import metadata
states <- read_sf("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/US_State_Boundaries/US_State_Boundaries.shp")
## remove non continental US states
states <- states %>% filter(NAME!="District of Columbia")
states <- states %>% filter(NAME!="U.S. Virgin Islands")
states <- states %>% filter(NAME!="Puerto Rico")

## monthly average temp from 2010 onward 
#avg_temp_IL <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/illinois_avgtemp.csv")
#avg_temp_IN <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/indiana_avgtemp.csv")
avg_temp_IA <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/iowa_avgtemp.csv")
#avg_temp_MI <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/michigan_avgtemp.csv")
#avg_temp_MN <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/minnesota_avgtemp.csv")
#avg_temp_NE <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/nebraska_avgtemp.csv")
#avg_temp_NC <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/northcarolina_avgtemp.csv")
#avg_temp_OH <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/ohio_avgtemp.csv")
#avg_temp_PA <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/pennsylvania_avgtemp.csv")
#avg_temp_KS <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/average temp_2010/kansas_avgtemp.csv")

## monthly total rainfall from 2010 onward 
#rain_IL <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/illinois_rain.csv")
#rain_IN <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/indiana_rain.csv")
rain_IA <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/iowa_rain.csv")
#rain_MI <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/michigan_rain.csv")
#rain_MN <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/minnesota_rain.csv")
#rain_NE <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/nebraska_rain.csv")
#rain_NC <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/northcarolina_rain.csv")
#rain_OH <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/ohio_rain.csv")
#rain_PA <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/pennsylvania_rain.csv")
#rain_KS <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/total precipitation_2010/kansas_rain.csv")

## clade assignments made separately using BV-BRC
tip_clades <- read.csv("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_v1/clade_assignment_updated.csv")
## combined maximum clade credibility (MCC) tree from TreeAnnotator
mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_v1/mcc_1990_v1_150k.trees")
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
tip_dates <- read.table("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_v1/aligned_HA_1990_dates.txt",
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
tip_meta <- merge(tip_meta, tip_clades, by = "sequence_name", all.x = TRUE)

#clean clade assignments
tip_meta <- tip_meta %>% mutate(clade = case_when(
                                clade == "Other-Human-2020" ~ "Other human",
                                clade == "Other-Human-1970" ~ "Other human",
                                clade == "Other-Human-1970-like" ~ "Other human",
                                clade == "Other-Human-2000" ~ "Other human",
                                clade == "Other-Human-2000-like" ~ "Other human",
                                clade == "Other-Human-2010" ~ "Other human",
                                clade == "Other-Human-2010-like" ~ "Other human",
                                TRUE ~ clade))
tip_meta <- tip_meta %>% mutate(clade = case_when(
                                clade == "1990.4-like" ~ "1990.4like",
                                clade == "1990.4" ~ "1990.4like",
                                clade == "1990.4.b1" ~ "1990.4.b1b2",
                                clade == "1990.4.b2" ~ "1990.4.b1b2",
                                clade == "" ~ "Missing",
                                clade == "2010.1-like" ~ "2010.1like",
                                clade == "2010.1" ~ "2010.1like",
                                is.na(clade) ~ "Missing",
                                TRUE ~ clade))
#table(tip_meta_clades$clade, useNA="always")

## once ESS values (excluding coalescent, join & prior if necessary, but ideally all values) are above 200, 
## export and combine .trees files using LogCombiner, then put combined tree into TreeAnnotator for the MCC tree

## part one B: MCC tree 

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
tree_clade <- ggtree(mcc_tree) %<+% tip_meta +
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
tree_clade_pruned <- ggtree(mcc_tree_pruned) %<+% tip_meta +
  geom_tippoint(aes(color = clade), alpha = 0.5, size = 3) +
  scale_color_paletteer_d("ggsci::default_ucscgb") +
  theme_tree2() + coord_cartesian(xlim = c(120, 175)) +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) +
  guides(color = guide_legend(ncol = 1))

#tree_clade_pruned
#ggsave("C:/Users/cgait/OneDrive/Desktop/clades_pruned.jpeg",width=20,height=25,units=c("cm"),tree_clade_pruned)


## part two: homogeneous Brownian diffusion ----

## template XML produced by BEAUTi (contains all 4589 taxa + sequences)
xml_path  <- "C:/Users/cgait/OneDrive/Desktop/1990_v2/1990_USflu_tree.xml"
out_dir   <- "C:/Users/cgait/OneDrive/Desktop/1990_v2/clade_xmls"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## MCMC settings (adjust as needed)
chain_length   <- 100000000   # 100 million
log_every      <- 10000       # sample every 10 000
save_every     <- 100000      # checkpoint every 100 000

## Parse template XML 
##    Build lookup tables: taxon_id → decimal date, taxon_id → sequence
doc <- read_xml(xml_path)

## extract all <taxon> nodes inside <taxa>
taxa_nodes <- xml_find_all(doc, ".//taxa/taxon")
taxon_ids  <- xml_attr(taxa_nodes, "id")
taxon_dates <- sapply(taxa_nodes, function(nd) {
  xml_attr(xml_find_first(nd, "date"), "value")
})
names(taxon_dates) <- taxon_ids

## extract all <sequence> nodes – each has a <taxon idref="..."/> + raw text
seq_nodes <- xml_find_all(doc, ".//alignment/sequence")
seq_ids   <- sapply(seq_nodes, function(nd) {
  xml_attr(xml_find_first(nd, "taxon"), "idref")
})
seq_seqs  <- sapply(seq_nodes, function(nd) {
  trimws(xml_text(nd))
})
names(seq_seqs) <- seq_ids
#cat("Parsed", length(taxon_ids), "taxa and", length(seq_ids), "sequences from template.\n")


## Lat/lon lookup from state swine centers
 locations <- data.frame(
   state = c("Iowa","Minnesota","Illinois","Indiana","Ohio","North Carolina",
             "Pennsylvania","Michigan","Wisconsin","Kentucky","Arkansas",
             "Oklahoma","Texas","Missouri","Kansas","Colorado","Wyoming",
             "Utah","Nebraska","South Dakota",
             ## ── states added to cover all tips ──
             "Alabama","Arizona","California","Florida","Georgia",
             "Louisiana","Maryland","Montana","New Mexico","New York",
             "North Dakota","Oregon","South Carolina","Tennessee",
             "Virginia","West Virginia"),
   lat = c(42.783, 44.105, 40.383, 40.843, 40.613, 35.236,
           40.475, 42.952, 43.631, 37.232, 34.300, 36.666,
           31.542, 40.074, 39.636, 39.021, 41.509, 38.172,
           41.452, 43.374,
           ## ── added states ──
           32.806,  34.048,  36.778,  27.995,  33.040,
           31.169,  39.046,  46.813,  34.840,  42.165,
           47.528,  44.572,  33.836,  35.518,
           37.769,  38.468),
   lon = c(-94.428, -95.262, -89.271, -86.314, -83.635, -78.193,
           -76.939, -85.810, -90.184, -87.148, -93.872, -99.812,
           -96.576, -92.102, -96.882, -102.390, -104.463, -113.663,
           -97.284, -97.486,
           ## ── added states ──
           -86.791, -111.094, -119.418, -81.760, -83.644,
           -91.867,  -76.641, -110.362, -106.249, -74.949,
           -99.784, -122.071,  -81.164, -86.580,
           -78.170,  -80.955))


## Build a named-vector lookup  state → c(lat, lon)
loc_lat <- setNames(locations$lat, locations$state)
loc_lon <- setNames(locations$lon, locations$state)
 

 ## Expand locations to cover ALL states that appear in tip_meta.
 ## These are approximate state centroids — update with your swine center
 ## coordinates once you have them for these states.
 extra_states <- data.frame(
   state = c("Alabama","Arizona","California","Florida","Georgia",
             "Louisiana","Maryland","Montana","New Mexico","New York",
             "North Dakota","Oregon","South Carolina","Tennessee",
             "Virginia","West Virginia"),
   lat   = c(32.806, 34.048, 36.778, 27.995, 33.040,
             31.169, 39.046, 46.813, 34.840, 42.165,
             47.528, 44.572, 33.836, 35.518,
             37.769, 38.468),
   lon   = c(-86.791, -111.094, -119.418, -81.760, -83.644,
             -91.867, -76.641, -110.362, -106.249, -74.949,
             -99.784, -122.071, -81.164, -86.580,
             -78.170, -80.955))
 
 ## only add states not already in the table
 extra_states <- extra_states[!extra_states$state %in% names(loc_lat), ]
 if (nrow(extra_states) > 0) {
   cat("Adding centroid coordinates for", nrow(extra_states), "extra states:",
       paste(extra_states$state, collapse = ", "), "\n")
   loc_lat <- c(loc_lat, setNames(extra_states$lat, extra_states$state))
   loc_lon <- c(loc_lon, setNames(extra_states$lon, extra_states$state))
 }
 
 ## ── 2b. Check which states in tip_meta are still un-geocoded ────────────────
 all_states_in_data <- unique(na.omit(tip_meta$state))
 missing_states     <- setdiff(all_states_in_data, names(loc_lat))
 if (length(missing_states) > 0) {
   warning("These states have no coordinates — their tips will be DROPPED ",
           "from clade XMLs:\n  ", paste(missing_states, collapse = ", "),
           "\n  → Add them to the locations data frame to include them.")
 }
 
 

## Generate one XML per clade
 ## filter tip_meta to only samples with a clade assignment
 tip_meta_assigned <- tip_meta %>% filter(!clade %in% c("Missing", NA_character_))
 ## get unique clades
 clades_to_run <- sort(unique(tip_meta_assigned$clade))
 #cat("Clades to process:", paste(clades_to_run, collapse = ", "), "\n\n")
 ## loop and write
 summary_rows <- list()
 
 for (cl in clades_to_run) {
   tips_in_clade <- tip_meta_assigned %>%
     filter(clade == cl) %>%
     pull(sequence_name)
   cat("── Clade:", cl, " (", length(tips_in_clade), " tips) ──\n")
   
   result <- build_clade_xml(
     clade_name   = cl,
     clade_tips   = tips_in_clade,
     tip_meta_df  = tip_meta,
     taxon_dates  = taxon_dates,
     seq_seqs     = seq_seqs,
     loc_lat      = loc_lat,
     loc_lon      = loc_lon,
     chain_length = chain_length,
     log_every    = log_every,
     save_every   = save_every)
   
   if (is.null(result)) next
   
   out_file <- file.path(out_dir, paste0(result$file_stem, ".xml"))
   writeLines(result$xml, out_file, useBytes = TRUE)
   cat("  → wrote", out_file, "\n")
   cat("    tips in XML:", result$n_tips,
       " | dropped (no coords):", result$n_dropped, "\n")
   
   summary_rows[[cl]] <- data.frame(
     clade            = cl,
     n_tips           = result$n_tips,
     n_dropped        = result$n_dropped,
     xml_file         = basename(out_file),
     stringsAsFactors = FALSE)
 }

## summary table
clade_summary <- bind_rows(summary_rows)
write.csv(clade_summary, file.path(out_dir, "clade_xml_summary.csv"), row.names = FALSE)

remove(doc, taxa_nodes, tree_clade, chain_length, log_every, out_dir, save_every, xml_path, tree_clade_pruned, 
       missing_clades, missing_in_tree, extra_states, cl, result)


## part two B: output from continuous lat/long diffusion ?? 



## part three: time-series outcomes ----

## create time-series prevalence data for each clade within each state
## add a year-month column for grouping
tip_meta_assigned <- tip_meta_assigned %>% mutate(year_month = floor_date(as.Date(date), "month"))

## overall monthly clade prevalence ---
overall_prev <- tip_meta_assigned %>% group_by(year_month, clade) %>% 
  summarise(count = n(), .groups = "drop") %>% group_by(year_month) %>%
  mutate(total = sum(count), prevalence = count / total) %>% ungroup()

## state-level monthly clade prevalence ---
state_prev <- tip_meta_assigned %>% group_by(year_month, state, clade) %>%
  summarise(count = n(), .groups = "drop") %>% group_by(year_month, state) %>% 
  mutate(total = sum(count), prevalence = count / total) %>% ungroup()

## define the clades and time window of interest
clades_of_interest <- c("2010.1like", "1990.4.a")
states_of_interest <- c("Iowa","Nebraska","Missouri","North Carolina")
date_start <- as.Date("2010-01-01")
date_end   <- as.Date("2026-12-31")

state_prev_filtered <- state_prev %>%
  filter(clade %in% clades_of_interest, year_month >= date_start, year_month <= date_end)
state_prev_filtered <- state_prev_filtered %>% filter(state %in% states_of_interest)

## plot clade prevalence over time in select states
state_clades <- ggplot(state_prev_filtered,
       aes(x = year_month, y = prevalence, color = clade)) +
geom_point(size = 4, alpha = 0.6) + scale_color_manual(values=c("pink3","orange2")) +
  facet_wrap(~ state, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Monthly Clade Prevalence by State (2010–2026)",
       x = "Month", y = "Prevalence (proportion of monthly samples)", color = "Clade") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
state_clades


## part four: maps ----

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

## swine centers for each state
## pulled just looking at the ESRI map lol will do more formally later
state_centers <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

US_samples <- ggplot() + geom_sf(data = states, aes(fill = n_sequences)) +
              geom_sf(data = state_centers, color="pink4", fill="gold", size=3, shape=22) +
              scale_fill_scico(palette = "hawaii", direction = -1) +
              theme_classic() 
#US_samples


## attach climate data over time to visualize? 

