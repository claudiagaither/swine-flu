## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!

## part one: data & baseline MCC tree ----

source("C:/Users/cgait/OneDrive/Desktop/swine flu/US_flu_functions.R")

## part one A: BEAST input
library(ape)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ggtree)
library(lubridate)
library(MASS)
library(msa)
library(paletteer)
library(Polychrome)
library(phangorn)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(scico)
#library(seraphim)
library(sf)
library(ShortRead)
library(stringr)
library(terra)
library(tidyverse)
library(treeio)
library(viridis)
library(writexl)
library(xml2)


## state-level metadata
states <- read_sf("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/US_State_Boundaries/US_State_Boundaries.shp")
## remove non continental US states
states <- states %>% filter(NAME!="District of Columbia")
states <- states %>% filter(NAME!="U.S. Virgin Islands")
states <- states %>% filter(NAME!="Puerto Rico")

## census of agriculture data
ag_census01 <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/agriculture data/ag_census/agcensus_alabama_hawaii.csv")


## climate data by state
base_path <- "C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/NOAA temp & rain/"
state_ids <- c(
  IL = "illinois",    IN = "indiana",       IA = "iowa",
  MI = "michigan",    MN = "minnesota",     NE = "nebraska",
  NC = "northcarolina", OH = "ohio",        PA = "pennsylvania",
  KS = "kansas",      AZ = "arizona",       AK = "arkansas",
  CO = "colorado",    KY = "kentucky",      MO = "missouri",
  MT = "montana",     NM = "newmexico",     OK = "oklahoma",
  SD = "southdakota", TX = "texas",         UT = "utah",
  WY = "wyoming")

## Import monthly average temperature across states
avg_temp <- lapply(state_ids, function(s) {
  read.csv(file.path(base_path, "average temp", paste0(s, "_avgtemp.csv")),
           skip = 2, header = TRUE)
})
names(avg_temp) <- names(state_ids)

## Import monthly total rainfall at state level ---
rain <- lapply(state_ids, function(s) {
  read.csv(file.path(base_path, "total precipitation", paste0(s, "_rain.csv")),
           skip = 2, header = TRUE)
})
names(rain) <- names(state_ids)

temp_df <- imap_dfr(avg_temp, function(df, abbr) {
  df$state <- abbr
  df$year_month <- as.Date(paste0(df$Date, "01"), format = "%Y%m%d")
  df})

rain_df <- imap_dfr(rain, function(df, abbr) {
  df$state <- abbr
  df$year_month <- as.Date(paste0(df$Date, "01"), format = "%Y%m%d")
  df})

## Filter to match the clade plot time range 
temp_df <- temp_df %>% filter(year_month >= as.Date("1997-01-01"))
rain_df <- rain_df %>% filter(year_month >= as.Date("1997-01-01"))


## clade assignments made separately using BV-BRC
tip_clades <- read.csv("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990 US flu/1990_v1/clade_assignment_updated.csv")
## combined maximum clade credibility (MCC) tree from TreeAnnotator
mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990 US flu/1990_v1/mcc_1990_v1_150k.trees")
options(ignore.negative.edge = TRUE)


## tip/taxa metadata
tip_dates <- read.table("C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990 US flu/1990_v1/aligned_HA_1990_dates.txt",
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
#length(intersect(tree_tips, clade_names))

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
#tree_clade_pruned <- ggtree(mcc_tree_pruned) %<+% tip_meta +
#  geom_tippoint(aes(color = clade), alpha = 0.5, size = 3) +
#  scale_color_paletteer_d("ggsci::default_ucscgb") +
#  theme_tree2() + coord_cartesian(xlim = c(120, 175)) +
#  theme(legend.position = "right", legend.text = element_text(size = 12),
#  legend.key.size = unit(0.8, "cm")) + guides(color = guide_legend(ncol = 1))
#tree_clade_pruned
#ggsave("C:/Users/cgait/OneDrive/Desktop/clades_pruned.jpeg",width=20,height=25,units=c("cm"),tree_clade_pruned)


## part two: random-walk diffusion xmls ----

## template XML produced by BEAUTi (contains all 4589 taxa + sequences)
#xml_path  <- "C:/Users/cgait/OneDrive/Desktop/1990_v2/1990_USflu_tree.xml"
#out_dir   <- "C:/Users/cgait/OneDrive/Desktop/clade_xmls"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## MCMC settings (adjust as needed)
#chain_length   <- 100000000   # 100 million
#log_every      <- 10000       # sample every 10 000
#save_every     <- 100000      # checkpoint every 100 000

## Parse template XML 
##    Build lookup tables: taxon_id → decimal date, taxon_id → sequence
#doc <- read_xml(xml_path)

## extract all <taxon> nodes inside <taxa>
#taxa_nodes <- xml_find_all(doc, ".//taxa/taxon")
#taxon_ids  <- xml_attr(taxa_nodes, "id")
#taxon_dates <- sapply(taxa_nodes, function(nd) {
#  xml_attr(xml_find_first(nd, "date"), "value")
#})
#names(taxon_dates) <- taxon_ids

## extract all <sequence> nodes – each has a <taxon idref="..."/> + raw text
#seq_nodes <- xml_find_all(doc, ".//alignment/sequence")
#seq_ids   <- sapply(seq_nodes, function(nd) {
#  xml_attr(xml_find_first(nd, "taxon"), "idref")
#})
#seq_seqs  <- sapply(seq_nodes, function(nd) {
#  trimws(xml_text(nd))
#})
#names(seq_seqs) <- seq_ids
#cat("Parsed", length(taxon_ids), "taxa and", length(seq_ids), "sequences from template.\n")


## Lat/lon lookup from state swine centers
 locations <- data.frame(
   state = c("Iowa","Minnesota","Illinois","Indiana","Ohio","North Carolina",
             "Pennsylvania","Michigan","Wisconsin","Kentucky","Arkansas",
             "Oklahoma","Texas","Missouri","Kansas","Colorado","Wyoming",
             "Utah","Nebraska","South Dakota",
             ## fewer than 10 sequences per state
             "Alabama","Arizona","California","Florida","Georgia",
             "Louisiana","Maryland","Montana","New Mexico","New York",
             "North Dakota","Oregon","South Carolina","Tennessee",
             "Virginia","West Virginia"),
   lat = c(42.783, 44.105, 40.383, 40.843, 40.613, 35.236,
           40.475, 42.952, 43.631, 37.232, 34.300, 36.666,
           31.542, 40.074, 39.636, 39.021, 41.509, 38.172,
           41.452, 43.374,
           32.806,  34.048,  36.778,  27.995,  33.040,
           31.169,  39.046,  46.813,  34.840,  42.165,
           47.528,  44.572,  33.836,  35.518,
           37.769,  38.468),
   lon = c(-94.428, -95.262, -89.271, -86.314, -83.635, -78.193,
           -76.939, -85.810, -90.184, -87.148, -93.872, -99.812,
           -96.576, -92.102, -96.882, -102.390, -104.463, -113.663,
           -97.284, -97.486,
           -86.791, -111.094, -119.418, -81.760, -83.644,
           -91.867,  -76.641, -110.362, -106.249, -74.949,
           -99.784, -122.071,  -81.164, -86.580,
           -78.170,  -80.955))
#write_xlsx(locations, "C:/Users/cgait/OneDrive/Desktop/locations.xlsx")

## Build a named-vector lookup  state → c(lat, lon)
loc_lat <- setNames(locations$lat, locations$state)
loc_lon <- setNames(locations$lon, locations$state)
 
 ## Check which states in tip_meta are still un-geocoded 
 #all_states_in_data <- unique(na.omit(tip_meta$state))
 #missing_states     <- setdiff(all_states_in_data, names(loc_lat))
 #if (length(missing_states) > 0) {
#   warning("These states have no coordinates — their tips will be DROPPED ",
#           "from clade XMLs:\n  ", paste(missing_states, collapse = ", "),
#           "\n  → Add them to the locations data frame to include them.")
# }
 
## Generate one XML per clade
## filter tip_meta to only samples with a clade assignment
tip_meta_assigned <- tip_meta %>% filter(!clade %in% c("Missing", NA_character_))
## get unique clades
#clades_to_run <- sort(unique(tip_meta_assigned$clade))
#cat("Clades to process:", paste(clades_to_run, collapse = ", "), "\n\n")
## loop and write
#summary_rows <- list()
#for (cl in clades_to_run) {
#   tips_in_clade <- tip_meta_assigned %>%
#     filter(clade == cl) %>%
#     pull(sequence_name)
#   cat("── Clade:", cl, " (", length(tips_in_clade), " tips) ──\n")
#   
#   result <- build_clade_xml(
#     clade_name   = cl,
#     clade_tips   = tips_in_clade,
#     tip_meta_df  = tip_meta,
#     taxon_dates  = taxon_dates,
#     seq_seqs     = seq_seqs,
#     loc_lat      = loc_lat,
#     loc_lon      = loc_lon,
#     chain_length = chain_length,
#     log_every    = log_every,
#     save_every   = save_every)
#  
#   if (is.null(result)) next
#   out_file <- file.path(out_dir, paste0(result$file_stem, ".xml"))
#   writeLines(result$xml, out_file, useBytes = TRUE)
#   cat("  → wrote", out_file, "\n")
#   cat("    tips in XML:", result$n_tips,
#       " | dropped (no coords):", result$n_dropped, "\n")
#   
#   summary_rows[[cl]] <- data.frame(
#     clade            = cl,
#     n_tips           = result$n_tips,
#     n_dropped        = result$n_dropped,
#     xml_file         = basename(out_file),
#     stringsAsFactors = FALSE)
# }

## summary table
#clade_summary <- bind_rows(summary_rows)
#write.csv(clade_summary, file.path(out_dir, "clade_xml_summary.csv"), row.names = FALSE)

remove(doc, taxa_nodes, tree_clade, chain_length, log_every, out_dir, save_every, xml_path, tree_clade_pruned, 
       missing_clades, missing_in_tree, extra_states, cl, result, rain, avg_temp)


## part two B: continuous-relaxed-random-walk surface



## part two C: branch-vector plot?
## cumulative Euclidean displacement vs phylogenetic distance from root/node age proxy



## part three: time-series data ----

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
clades_of_interest <- c("2010.1like", "1990.4.a","1990.1")
states_of_interest <- c("Iowa","Nebraska","Missouri","North Carolina")
date_start <- as.Date("1990-01-01")
date_end   <- as.Date("2026-4-01")
state_prev_filtered <- state_prev %>% filter(state %in% states_of_interest)
#state_prev_filtered <- state_prev_filtered %>% filter(clade %in% clades_of_interest)

state_lookup <- c(IA = "Iowa", NE = "Nebraska", MO = "Missouri", NC = "North Carolina")
temp_df$state <- state_lookup[temp_df$state]
rain_df$state <- state_lookup[rain_df$state]
temp_df <- temp_df %>% filter(!is.na(state))
rain_df <- rain_df %>% filter(!is.na(state))

## plot clade prevalence over time in select states
## overlay monthly rainfall for these states over the same time?? as a line
# Choose a scaling factor to map temperature onto the prevalence axis

rain_scale <- max(state_prev_filtered$prevalence, na.rm = TRUE) / 
  max(rain_df$Value, na.rm = TRUE)
state_clades_rain <- ggplot(state_prev_filtered,
  aes(x = year_month, y = prevalence, color = clade)) +
  geom_point(size = 3.5, alpha = 0.6) +
  geom_line(data = rain_df, aes(x = year_month, y = Value * rain_scale, color = NULL),
            color = "steelblue", linewidth = 0.5, alpha = 0.5) +
#  scale_color_manual(values = c("skyblue", "orange2", "pink")) +
  facet_wrap(~ state, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format(),
    sec.axis = sec_axis(~ . / rain_scale, name = "Total Rainfall (in)")) +
  labs(title = "Monthly Clade Prevalence & Rainfall by State (1997–2026)",
       x = "Month", y = "Prevalence", color = "Clade") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
#state_clades_rain

## outcome of a new dominant clade over time for each state
dominant <- state_prev_filtered %>% group_by(state, year_month) %>%
  slice_max(prevalence, n = 1, with_ties = FALSE) %>%
  ungroup() %>% arrange(state, year_month) %>% group_by(state) %>%
  mutate(prev_dominant = lag(clade),new_dominant = ifelse(
    clade != prev_dominant & !is.na(prev_dominant), 1, 0)) %>% ungroup()


## code indicator for major clades
dominant <- dominant %>% mutate(major_clade = case_when(
            clade == "1990.1" ~ "1990", clade == "1990.4.a" ~ "1990",
            clade == "1990.4.b1b2" ~ "1990", clade == "1990.4.d" ~ "1990",
            clade == "1990.4.f" ~ "1990", clade == "1990.4.i" ~ "1990",
            clade == "1990.4.k" ~ "1990", clade == "1990.4like" ~ "1990",
            clade == "2010.1like" ~ "2010", clade == "Other human" ~ "Human",
            TRUE ~ NA))

## indicator outcome of shift between the 2 major clades 
dominant <- dominant %>% arrange(state, year_month) %>%
  group_by(state) %>% mutate(prev_major = lag(major_clade),
    major_dominant = ifelse(major_clade != prev_major &
        !is.na(prev_major) & major_clade %in% c("1990", "2010") &
    prev_major %in% c("1990", "2010"), 1, 0)) %>% ungroup()

shift_events <- dominant %>% filter(major_dominant == 1)
#shift_events <- dominant %>% filter(new_dominant == 1)

state_shifts_rain <- ggplot(state_prev_filtered, aes(x = year_month, y = prevalence, color = clade)) +
  geom_line(data = rain_df, aes(x = year_month, y = Value * rain_scale, color = NULL),
            color = "cornflowerblue", linewidth = 0.7, alpha = 0.6) +
  geom_vline(data = shift_events, aes(xintercept = year_month), color = "pink3", alpha = 0.7) +
  facet_wrap(~ state, scales = "free_y") + scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format(),
                     sec.axis = sec_axis(~ . / rain_scale, name = "Total Rainfall (in)")) +
  labs(title = "Clade Replacement Events & Rainfall by State", x="Month", y="Prevalence", color="Clade") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
#state_shifts_rain



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
#states <- states %>% filter(n_sequences>=10)

## swine centers for each state
## pulled just looking at the ESRI map lol will do more formally later
state_centers <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

US_samples <- ggplot() + geom_sf(data = states, aes(fill = n_sequences)) +
              geom_sf(data = state_centers, color="pink4", fill="gold", size=3, shape=22) +
              scale_fill_scico(palette = "hawaii", direction = -1) +
              theme_classic() 
#US_samples


## KDE for node density surfaces within subclades??
tree_files <- list(
  "1990.4.a"   = "C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_1990.4.a_v2/mcc_1990.4.a_v2.trees",
  "2010.1like" = "C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_2010.1like_v2/mcc_1990_2010.1like.trees",
  "2010.2"     = "C:/Users/cgait/OneDrive/Desktop/2010.2_v1/2010.2_v1.trees",
  "1990.4.b"   = "C:/Users/cgait/OneDrive/Desktop/1990.4b_v1/1990.4.b1b2_v1.trees")


## Returns a named list of ggplot objects, and writes a PNG for each.
#diffusion_maps <- mapply(
#  FUN = function(file, name) {
#    make_diffusion_map(
#      tree_file     = file,
#      subclade_name = name,
#      save_path     = paste0(name, "_surface.png"))
#  },
#  file = tree_files,
#  name = names(tree_files),
#  SIMPLIFY = FALSE)

## Access individually:
#diffusion_maps[["1990.4.a"]]
#diffusion_maps[["2010.1like"]]
#make_diffusion_map(tree_files[["2010.2"]], "2010.2", 
#                   bw_mult_lon = 1.0, bw_mult_lat = 3.0)
#diffusion_maps[["1990.4.b"]]

#ggsave("2010.1like_surface.png", plot=diffusion_2010.1like, width=12, height=8, dpi=300, bg="white")
