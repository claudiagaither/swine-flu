## H3N2 influenza A transmission among domestic swine in the US
## phylogeography and ecological predictors of transmission!

## part one:    metadata ----
set.seed(1738)
source("C:/Users/cgait/OneDrive/Desktop/swine flu/US_flu_functions.R")
library(ape)
library(Biostrings)
library(climateR)
library(conflicted)
library(data.table)
library(dplyr)
library(exactextractr)
library(forecast) 
library(geepack)
library(ggnewscale)
library(ggplot2)
library(ggthemes)
library(ggtree)
library(gt)
library(lubridate)
library(MASS)
library(msa)
library(paletteer)
library(patchwork)
library(Polychrome)
library(phangorn)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(scico)
library(sf)
library(ShortRead)
library(stringr)
library(terra)
library(tidyr)
library(tidyverse)
library(tigris)
library(tseries)  
library(treeio)
library(trend)
library(viridis)
library(writexl)
library(xml2)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count",  "dplyr")
conflict_prefer("summarise", "dplyr")
conflicts_prefer(base::as.data.frame)
conflicts_prefer(base::unique)


## import & clean metadata
states <- read_sf("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/US_State_Boundaries/US_State_Boundaries.shp")
hog_census <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/agriculture data/hog_census.csv")
hog_survey_01 <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/agriculture data/ag_survey_alabama_mississippi.csv")
hog_survey_02 <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/agriculture data/ag_survey_missouri_tennessee.csv")
hog_survey_03 <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/agriculture data/ag_survey_texas_wisconsin.csv")
climate_state_wt <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate data/climate_state_wt.csv")

## remove non continental US states
states <- states %>% filter(NAME!="District of Columbia")
states <- states %>% filter(NAME!="U.S. Virgin Islands")
states <- states %>% filter(NAME!="Puerto Rico")

## filter for counties represented by hog census?
hog_survey <- bind_rows(hog_survey_01, hog_survey_02, hog_survey_03)
county_names <-hog_survey[,c("State","County","Value")]
hogs_county <- hog_census[,c("STATE_ANSI","COUNTY_ANSI","STATE_NAME")]

## per-county hog INVENTORY weights from the annual survey
survey_weights <- hog_survey %>% filter(Geo.Level  == "COUNTY",
                                        Commodity  == "HOGS",
         Data.Item == "HOGS - INVENTORY",   # inventory, not sales
         Domain == "TOTAL", Domain.Category == "NOT SPECIFIED",
         !is.na(State.ANSI), !is.na(County.ANSI)) %>%   # drops "OTHER COUNTIES"/district rows
  transmute(
    GEOID       = sprintf("%02d%03d", as.integer(State.ANSI), as.integer(County.ANSI)),
    survey_year = as.integer(Year),
    head        = as.numeric(gsub(",", "", Value))      # "(D)"/"(Z)" -> NA on coercion
  ) %>% filter(!is.na(head)) %>% group_by(GEOID, survey_year) %>%     # collapse any duplicate Periods in a year
  summarise(head = mean(head, na.rm = TRUE), .groups = "drop")

## suppression check: usable counties per year
survey_weights %>% dplyr::count(survey_year) %>% arrange(survey_year)


## swine center selection algorithm based on Census of Agriculture data
## Get county centroids for all US counties
counties_sf <- tigris::counties(cb = TRUE, progress_bar = FALSE)

## Extract lat/lon from county geometries
county_coords <- counties_sf %>% st_centroid() %>% mutate(GEOID = GEOID,
    county_lat = st_coordinates(.)[, "Y"],
    county_lon = st_coordinates(.)[, "X"]) %>%
  st_drop_geometry() %>% select(GEOID, county_lat, county_lon)

## Average hog weights across all survey years (1998-2025) by county
county_hog_avg <- survey_weights %>%
  group_by(GEOID) %>% summarise(avg_hog_count = mean(head, na.rm = TRUE), 
            .groups = "drop")

## Join coordinates with hog weights
county_data <- county_coords %>% inner_join(county_hog_avg, by = "GEOID")

## Get state FIPS codes to match your GEOID format
state_fips <- tigris::states(cb = TRUE) %>%
  st_drop_geometry() %>% select(STATEFP, NAME) %>%
  dplyr::rename(state = NAME) %>% mutate(state_fips = as.numeric(STATEFP))

## need lat and lon columns
## there are also sequences from Florida, Louisiana, Maryland, Utah & Wyoming
## but no survey data apparently?
locations <- data.frame(state = c("Iowa","Minnesota","Illinois","Indiana","Ohio","North Carolina",
                        "Pennsylvania","Michigan","Wisconsin","Kentucky","Arkansas",
                        "Oklahoma","Texas","Missouri","Kansas","Colorado","Nebraska","South Dakota",
                        # fewer than 10 sequences per state
                        "Alabama","Arizona","California","Georgia","Montana","New Mexico","New York",
                        "North Dakota","Oregon","South Carolina","Tennessee","Virginia","West Virginia"))

## Calculate weighted-average coordinates by state
state_locations <- county_data %>% mutate(state_fips = as.numeric(substr(GEOID, 1, 2))) %>%
  left_join(state_fips, by = "state_fips") %>% filter(state %in% locations$state) %>%  # filter to your states of interest
  group_by(state) %>% summarise(
    latitude = weighted.mean(county_lat, avg_hog_count, na.rm = TRUE),
    longitude = weighted.mean(county_lon, avg_hog_count, na.rm = TRUE),
    total_hogs = sum(avg_hog_count, na.rm = TRUE),
    n_counties = n(), .groups = "drop") %>% arrange(state)

## Convert to named vectors for lookup (as you were planning)
loc_lat <- setNames(state_locations$latitude, state_locations$state)
loc_lon <- setNames(state_locations$longitude, state_locations$state)


## part two:    MCC tree and tips ----

## combined maximum clade credibility (MCC) tree from TreeAnnotator
mcc_tree <- read.beast("C:/Users/cgait/OneDrive/Desktop/swine flu/sequence data/mcc_1990_v3.trees")
options(ignore.negative.edge = TRUE)

## tip/taxa metadata (dates, states & BV_BRC clade assignments)
tip_meta <- read.table("C:/Users/cgait/OneDrive/Desktop/swine flu/sequence data/H3_1990-2026_metadata_final.tsv",
                        header = TRUE, sep = "\t")
tip_meta$date <- as.Date(tip_meta$date)

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
#table(tip_meta$clade, useNA="always")

## once ESS values (excluding coalescent, join & prior if necessary, but ideally all values) are above 200, 
## export and combine .trees files using LogCombiner, then put combined tree into TreeAnnotator for the MCC tree

# extract and view data
tree_phylo <- mcc_tree@phylo

#color tips by clade assignments made using BV-BRC
#tree_clade <- ggtree(mcc_tree) %<+% tip_meta +
#  geom_tippoint(aes(color = clade), alpha=0.5, size = 3) +
#  scale_color_paletteer_d("ggsci::default_ucscgb") + 
#  theme_tree2() + coord_cartesian(xlim = c(1990, 2033)) +
#  theme(legend.position = "right", legend.text = element_text(size = 12),
#        legend.key.size = unit(0.8, "cm")) + guides(color = guide_legend(ncol = 1))
#tree_clade

# See how many sequences match between tree and clade assignments
tree_tips <- mcc_tree@phylo$tip.label
clade_names <- tip_meta$sequence_name
length(base::intersect(tree_tips, clade_names))

# See which clade assignments are missing (not matching tree tips)
missing_in_tree <- base::setdiff(clade_names, tree_tips)

# See which tree tips don't have clade assignments
missing_clades <- base::setdiff(tree_tips, clade_names)

norm_key <- function(x) {
  x <- gsub("_dup[0-9]+$", "", x)   # drop dedup suffixes
  x <- gsub(",", "", x)             # drop commas the tree software removed
  x <- gsub("\\s+", " ", x)         # collapse whitespace
  trimws(x)
}

lookup <- tip_meta %>% mutate(key = norm_key(sequence_name)) %>% distinct(key, .keep_all = TRUE) %>%      # handles the 17 dup rows
  select(key, clade)

tree_df <- data.frame(tip = tree_tips, key = norm_key(tree_tips)) %>% left_join(lookup, by = "key")
#sum(is.na(tree_df$clade))     # residual unmatched tips, n=50

missing_clades <- tree_df$tip[is.na(tree_df$clade)]
#length(missing_clades) 

## prune missing_clades tips from tree and re-plot?
mcc_tree_pruned <- treeio::drop.tip(mcc_tree, missing_clades)

# sanity check: should drop by exactly 50
#ape::Ntip(mcc_tree@phylo)            # 4589
#ape::Ntip(mcc_tree_pruned@phylo)     # 4539

# Re-plot with pruned tree
clade_df <- tree_df[, c("tip", "clade")]
#clades_pruned <- ggtree(mcc_tree_pruned) %<+% clade_df + geom_tippoint(aes(color = clade))
#clades_pruned
#ggsave("C:/Users/cgait/OneDrive/Desktop/clades_pruned.jpeg",width=20,height=25,units=c("cm"), clades_pruned)


## part three:  random-walk diffusion xmls ----

## template XML produced by BEAUTi (contains all 4589 taxa + sequences)
#xml_path  <- "C:/Users/cgait/OneDrive/Desktop/BEAST runs/1990_v3/1990_USflu_tree.xml"
#out_dir   <- "C:/Users/cgait/OneDrive/Desktop/clade_xmls"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## need to filter full mcc template to remove states with < 30 sequences
## then randomly sample the minimum number (should be 31) of sequences from each state
## and THEN make the clade-specific xmls from this downsampled set of sequences

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

## Check which states in tip_meta are still un-geocoded 
#all_states_in_data <- unique(na.omit(tip_meta$state))
#missing_states     <- base::setdiff(all_states_in_data, names(loc_lat))
#if (length(missing_states) > 0) {
#   warning("These states have no coordinates — their tips will be DROPPED ",
#         "from clade XMLs:\n  ", paste(missing_states, collapse = ", "),
#           "\n  → Add them to the locations data frame to include them.")
# }
 
## Generate one XML per clade
## filter tip_meta to only samples with a clade assignment
tip_meta_assigned <- tip_meta %>% filter(!clade %in% c("unassigned", NA_character_))
tip_meta_assigned <- tip_meta_assigned %>% distinct(sequence_name, .keep_all = TRUE)

## get unique clades
#clades_to_run <- sort(unique(tip_meta_assigned$clade))
#cat("Clades to process:", paste(clades_to_run, collapse = ", "), "\n\n")

## Downsample sequences: each state contributes equal n of sequences 
## Apply the SAME filters build_clade_xml uses internally so the per-state
## minimum reflects sequences that will actually make it into an XML.
#pool <- tip_meta_assigned %>% filter(sequence_name %in% names(taxon_dates),
#         sequence_name %in% names(seq_seqs), !is.na(state), state %in% names(loc_lat))
## Drop states with < 30 sequences across the whole (filtered) pool
#min_state_n  <- 30
#state_counts <- pool %>% dplyr::count(state, name = "n")
#keep_states     <- state_counts$state[state_counts$n >= min_state_n]
#dropped_states  <- state_counts$state[state_counts$n <  min_state_n]
#cat("State counts in pool (sorted):\n"); print(state_counts %>% arrange(n))
#cat("\nDropping", length(dropped_states), "states with <", min_state_n,
#   "seqs:\n  ", paste(dropped_states, collapse = ", "), "\n")
#pool <- pool %>% filter(state %in% keep_states)
### Take the minimum surviving state count -> per-state sample size
#n_per_state <- min(state_counts$n[state_counts$state %in% keep_states])
#cat("\nDownsampling each kept state to n =", n_per_state, "sequences.\n")
## Random sample, equal n from every kept state
#downsampled <- pool %>% group_by(state) %>% slice_sample(n = n_per_state) %>% ungroup()
#allowed_tips <- downsampled$sequence_name
#cat("Total seqs after downsampling:", length(allowed_tips),
#    "(", length(keep_states), "states x", n_per_state, "per state)\n\n")
## sanity check
#downsampled_clades <- downsampled %>% dplyr::count(clade, state) # per-clade breakdown

## loop and write xmls
#summary_rows <- list()
#for (cl in clades_to_run) {
#  tips_in_clade <- tip_meta_assigned %>%
#    filter(clade == cl, sequence_name %in% allowed_tips) %>%
#    pull(sequence_name)
#  cat("── Clade:", cl, " (", length(tips_in_clade), " tips) ──\n")
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
#write.csv(downsampled %>% select(sequence_name, state, clade, date),
#          file.path(out_dir, "downsampled_taxa.csv"), row.names = FALSE)

remove(doc, taxa_nodes, chain_length, log_every, out_dir, save_every, xml_path, 
       missing_clades, missing_in_tree, extra_states, cl, result, rain, avg_temp, hog_survey_01, hog_survey_02,
       hog_survey_03, mcc_tree, tip_dates, tree_phylo, hog_census)


## part four:   ARIMA(X)? ----

## create time-series prevalence data for each clade within each state
## add a year-month column for grouping
tip_meta_assigned <- tip_meta_assigned %>% mutate(year_month = floor_date(as.Date(date), "month"))
tip_meta_assigned$year <- as.integer(format(tip_meta_assigned$year_month, "%Y"))
tip_meta_assigned <- tip_meta_assigned %>% filter(year >= 2003)
#hist(tip_meta_assigned$year)
#table(tip_meta_assigned$state)

## overall monthly clade prevalence 
overall_prev <- tip_meta_assigned %>% group_by(year_month, clade) %>% 
  summarise(count = n(), .groups = "drop") %>% group_by(year_month) %>%
  mutate(total = sum(count), prevalence = count / total) %>% ungroup()

## state-level monthly clade prevalence 
state_prev <- tip_meta_assigned %>% group_by(year_month, state, clade) %>%
  summarise(count = n(), .groups = "drop") %>% group_by(year_month, state) %>% 
  mutate(total = sum(count), prevalence = count / total) %>% ungroup()

## define the clades and time window of interest
clades_of_interest <- c("2010.1like", "1990.4.a","1990.1")
states_of_interest <- c("Iowa","Nebraska","Missouri","North Carolina")
date_start <- as.Date("1990-01-01")
date_end   <- as.Date("2024-12-31")
state_prev_filtered <- state_prev %>% filter(state %in% states_of_interest)

## outcome of a new dominant clade over time for each state
dominant <- state_prev_filtered %>% group_by(state, year_month) %>%
  slice_max(prevalence, n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(state, year_month) %>% group_by(state) %>%
  mutate(prev_dominant = dplyr::lag(clade),new_dominant = ifelse(clade != prev_dominant & !is.na(prev_dominant), 1, 0)) %>% ungroup()

## code indicator for lineages
dominant <- dominant %>% mutate(major_clade = case_when(
            clade == "1990.1" ~ "1990", clade == "1990.4.a" ~ "1990",
            clade == "1990.4.b1b2" ~ "1990", clade == "1990.4.d" ~ "1990",
            clade == "1990.4.f" ~ "1990", clade == "1990.4.i" ~ "1990",
            clade == "1990.4.k" ~ "1990", clade == "1990.4like" ~ "1990",
            clade == "2010.1like" ~ "2010", clade == "Other human" ~ "Human", TRUE ~ NA))

## indicator outcome of shift between the 2 major lineages
dominant <- dominant %>% arrange(state, year_month) %>% group_by(state) %>% mutate(prev_major = dplyr::lag(major_clade),
    major_dominant = ifelse(major_clade != prev_major & !is.na(prev_major) & major_clade %in% c("1990", "2010") &
    prev_major %in% c("1990", "2010"), 1, 0)) %>% ungroup()

## change bins for precipitation and temp variables for forest plots
climate_state_wt$Precipitation <- climate_state_wt$ppt/50
climate_state_wt$Min_temp <- climate_state_wt$tmin/10
climate_state_wt$Max_temp <- climate_state_wt$tmax/10
climate_state_wt$Average_temp <- climate_state_wt$tavg/10
climate_state_wt$Vapor_pressure <- climate_state_wt$vap
climate_state_wt$Absolute_humidity <- climate_state_wt$abs_humidity 
climate_state_wt$Wind_speed <- climate_state_wt$ws

clim_vars <- c("Max_temp", "Min_temp", "Average_temp", "Precipitation", "Vapor_pressure", "Wind_speed", "Absolute_humidity")
#is_major_shift <- function(df) df$major_dominant == 1

## long-format climate data
clim <- climate_state_wt %>% mutate(date = as.Date(sprintf("%d-%02d-01", year, month)))
clim_long <- clim %>% pivot_longer(all_of(clim_vars), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = clim_vars))  # keeps facet order
clim_long <- clim_long %>% filter(date <= "2024-12-31")
clim_long <- clim_long %>% filter(date >= "2012-01-01")

## based on the trend summary, the slope for max temp and avg temp have p values under 0.05
## min temp has a p-value of 0.0544, but others are all over 0.10

# part three B: ARIMA parameterization sensitivity 
# Picks the differencing order d (and seasonal D) PROGRAMMATICALLY via unit-root
# tests, flags the choice, then sweeps only (p, q) at that fixed d/D. Because
# d and D are constant across every model, the AICc / ΔAICc column in the final
# output is fully comparable -- no cross-d apples-to-oranges problem.

ds <- as.Date("2012-01-01"); de <- as.Date("2024-12-01")
#run_2010 <- run_clade_arimax(state_prev, climate_state_wt, "Iowa", "2010.1like", ds, de)

#run_2010$predictor_results      # which climate var best explains the clade
#render_predictor_table(run_2010)
#render_sensitivity_table(run_2010)
#render_sensitivity_heatmap(run_2010)
#check_clade_model(run_2010)      # residual diagnostics on the winning order


## expand to multiple states and clades
combos <- expand_grid(state = c("Iowa","Illinois","Indiana","Ohio","Pennsylvania","North Carolina",
                      "Minnesota","South Dakota","Nebraska",
                      "Oklahoma","Missouri","Kansas","Arkansas"),
                      clade = c("2010.1like","1990.4.a"))
#runs <- pmap(combos, ~ run_clade_arimax(state_prev, climate_state_wt, ..1, ..2, ds, de))
#all_predictor_aic <- map_dfr(runs, "predictor_results")   # stacked, has state+clade cols

#predictor_summary <- build_predictor_summary(runs, adjust = TRUE)
#predictor_table <- render_significance_table(predictor_summary)
#write.csv(predictor_summary, "C:/Users/cgait/OneDrive/Desktop/summary_predictor.csv")
#summary_stability <- stability_summary(runs)
#write.csv(summary_stability, "C:/Users/cgait/OneDrive/Desktop/summary_stability.csv")

## based on the predictor summary for first four states (Iowa, Missouri, Nebraska, North Carolina)
## we see the onnly significant predictors in Missouri, but stable across parameterizations
## for minimum temp (derivative) and precpitation... n = 75 obs so not the largest even

#predictor_summary <- read.csv("C:/Users/cgait/OneDrive/Desktop/swine flu/climate models/summary_predictor.csv")
# geographic: Midwest swine belt grouped, North Carolina (Southeast) last
#forest_predictor(predictor_summary, state_levels = c("Arkansas","Illinois","Indiana","Iowa",
#                                   "Kansas", "Minnesota", "Missouri", "Nebraska", "North Carolina",                  
#                                   "Ohio", "Oklahoma", "Pennsylvania","South Dakota"))
#ggsave("predictor_forest.png", width = 14, height = 10, dpi = 300, bg = "white")


## part five:   gee for dominant clades ----

## Same climate predictors as part four (predictor_columns(): detrended temps
## *_dt + raw ppt/vap/ws/abs_humidity), but BINARY switching outcomes modeled
## by GEE (binomial / logit) to handle within-state temporal correlation:
##   new_dominant   : dominant clade changed vs the previous sequenced month
##   major_dominant : dominant major lineage switched between 1990 and 2010
## One model per state -> "top predictor per state" forest plot, faceted by
## outcome. corstr = exchangeable, robust (sandwich) SE; months clustered by
## calendar year. See part four-b of the functions file for the engine.


## window. NOTE: tip_meta_assigned is filtered to year >= 2003 upstream (part
## four), so although we ask for 2002 the series effectively begins 2003. To
## truly start in 2002, relax that filter where state_prev is built.
gee_start <- as.Date("2003-01-01")
gee_end   <- as.Date("2024-12-31")

## states with >= 20 sequences in the window
gee_states_tbl <- qualifying_states(state_prev, gee_start, gee_end, min_seqs = 30)
gee_states     <- gee_states_tbl$state
#print(gee_states_tbl)

## monthly dominant-clade panel for all qualifying states
## (both outcomes + climate predictors, temps detrended within state)
#dom_panel <- build_dominant_panel(state_prev, climate_state_wt, gee_states,
#                                  gee_start, gee_end)

## per-state, per-outcome GEE: rank climate predictors by QIC, keep the best,
## robust-Wald inference + odds ratios, BH-adjusted across the family
#gee_summary <- build_gee_summary(dom_panel, gee_states,
#                                 outcomes   = c("new_dominant", "major_dominant"),
#                                 id_var     = "year",
#                                 corstr     = "exchangeable",
#                                 min_events = 3,
#                                 adjust     = TRUE)
#write.csv(gee_summary, "C:/Users/cgait/OneDrive/Desktop/gee_summary.csv", row.names = FALSE)

## gt table + forest plot (one row per state, faceted by outcome)
#render_gee_table(gee_summary)
#forest_gee(gee_summary)
#ggsave("gee_predictor_forest.png", width = 15, height = 10, dpi = 300, bg = "white")


## optional cross-check: one properly month-ordered AR-1 GEE clustered on STATE
## per (outcome, predictor), to sanity-check the per-state estimates
#fit_gee_pooled(dom_panel, "new_dominant",   "tavg_dt", corstr = "ar1")
#fit_gee_pooled(dom_panel, "major_dominant", "ppt",     corstr = "ar1")


## part six:    maps ----

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

## filter states for inclusion
states <- states %>% filter(NAME %in% state_locations$state)
states <- states %>% filter(n_sequences>30)
state_locations <- state_locations %>% filter(state %in% states$NAME)

## swine centers for each state
state_centers <- st_as_sf(state_locations, coords = c("longitude", "latitude"), crs = 4326)

## KDE for node density surfaces within clades
tree_files <- list(
  "1990.4.a"= "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/1990.4.a/mcc_1990.4.a_downsampled.trees",
  "1990.4.b"= "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/1990.4.b1b2/mcc_1990.4.b1b2_downsampled.trees",
  "2010.1"  = "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/2010.1like/mcc_2010.1like_downsampled.trees",
  "2010.2"  = "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/2010.2/mcc_2010.2_downsampled.trees")
## smaller subclades (make tips = 13 have much lower densities so not included in results)
#  "1990.4.e"= "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/smaller_subclades/mcc_1990.4.e.trees",
#  "1990.4.f"= "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/smaller_subclades/mcc_1990.4.f.trees")
#  "1990.4.i"= "C:/Users/cgait/OneDrive/Desktop/BEAST runs/downsampled_subclades/smaller_subclades/mcc_1990.4.i.trees")

## compute all four KDE surfaces (trees read once each)
kde_list <- mapply(
  FUN           = compute_diffusion_kde,
  tree_file     = tree_files,
  subclade_name = names(tree_files),
  SIMPLIFY      = FALSE)

## global density range across all four
global_limits <- range(unlist(lapply(kde_list, `[[`, "density")))

## build maps on the shared scale
diffusion_maps <- mapply(
  FUN  = function(df, name) plot_diffusion_map(df, name, fill_limits = global_limits),
  df   = kde_list,
  name = names(kde_list),
  SIMPLIFY = FALSE)
#subclade_maps <- wrap_plots(diffusion_maps, nrow = 2, ncol = 2 + plot_layout(guides = "collect"))
#ggsave("C:/Users/cgait/OneDrive/Desktop/subclade_maps_smaller.jpeg", width=25,height=15,units=c("cm"), subclade_maps)

## node posterior support within subclades
posterior_list <- mapply(
  FUN           = compute_node_posteriors,
  tree_file     = tree_files,
  subclade_name = names(tree_files),
  SIMPLIFY      = FALSE)

# Proportion above a single threshold for every clade
sapply(posterior_list, function(df) mean(df$posterior > 0.5))
sapply(posterior_list, function(df) mean(df$posterior > 0.8))

## posterior is a probability -> fixed 0-1 scale across all four
posterior_maps <- mapply(
  FUN  = function(df, name) plot_posterior_map(df, name, fill_limits = c(0, 1)),
  df   = posterior_list,
  name = names(posterior_list),
  SIMPLIFY = FALSE)
#posterior_maps <- wrap_plots(posterior_maps, nrow = 2, ncol = 2 + plot_layout(guides = "collect"))
#ggsave("C:/Users/cgait/OneDrive/Desktop/posterior_maps_smaller.jpeg", width=25, height=15, units = c("cm"), posterior_maps)

