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
library(phangorn)
library(readxl)
library(ShortRead)
library(stringr)
library(tidyverse)
library(treeio)
library(writexl)


## part one: BEAST input ----

## remove problematic sequence: A/swine/Illinois/A00857131/2011|EPI_ISL_121898|A_/_H3N2||||2011-09-24|HA|4 from fasta before alignment
## this is what we will prune from the version 6 tree visualized below!

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
                                state == "Dakota_Superior" ~ "North Dakota",   ## let's double check this!
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
#log <- read.table("C:/Users/cgait/OneDrive/Desktop/BEAST runs/H3N2_2010_v6/comb_epiflu_HA2010_v6.log",
#                  header = TRUE, comment.char = "#", sep = "\t")

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
#tree_states <- ggtree(mcc_tree) %<+% tip_meta +
#  geom_tippoint(aes(color = state), size = 2) +
#  theme_tree2() + coord_cartesian(xlim = c(105, 135)) +
#  theme(legend.position = "right", legend.text = element_text(size = 8),
#        legend.key.size = unit(0.4, "cm")) + guides(color = guide_legend(ncol = 2))
#tree_states

#color tips by broader region
#add back in tip labels to export for qualitative clade assignment!
tree_region <- ggtree(mcc_tree) %<+% tip_meta +
  geom_tippoint(aes(color = region), alpha=0.5, size = 3) +
  scale_color_brewer(palette="Accent") + theme_tree2() + coord_cartesian(xlim = c(95, 135)) +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) + guides(color = guide_legend(ncol = 1))
#tree_region


## prune 2 outlier sequences from version 6 tree
remove_tips <- c("A/swine/Illinois/A00857131/2011|EPI_ISL_121898|A_/_H3N2||||2011-09-24|HA|4")

# verify both tips actually exist in the tree before pruning
missing <- setdiff(remove_tips, mcc_tree@phylo$tip.label)
if (length(missing)) stop(paste("These tips are not in the tree:", paste(missing, collapse = ", ")))

# prune tree — drop.tip removes the named tips directly
mcc_tree <- treeio::drop.tip(mcc_tree, remove_tips)

# remove dates that are not in the tree
tip_meta <- tip_meta[tip_meta$label %in% mcc_tree@phylo$tip.label, , drop = FALSE]
# verify
#stopifnot(length(mcc_tree@phylo$tip.label) == nrow(tip_meta))

# prune metadata to match
tip_meta <- tip_meta[!tip_meta$label %in% remove_tips, , drop = FALSE]

# verify tip count matches
#stopifnot(length(mcc_tree@phylo$tip.label) == nrow(tip_meta))


## export excel sheet of all taxa in the top to bottom order they appear on the tree
## if we want to manually edit for whatever reason
# Get tip order as displayed in ggtree (top to bottom)
#tip_order <- ggtree(mcc_tree)$data %>%
#  filter(isTip) %>% arrange(y) %>% pull(label)      # y-axis order = visual order
# Build a data frame for export
#clade_sheet <- data.frame(tree_order = seq_along(tip_order),
#  label = tip_order, clade = NA_character_) %>%       # you'll fill this in manually
#  left_join(tip_meta, by = "label")
#write_xlsx(clade_sheet, "manual_clade_assignments.xlsx")



## part three: baseline clades ----

## algorithmic assignment based on posterior probability support for nodes!
## clade assignment based on posterior probability and target clade size, steps 1-9

tree <- mcc_tree@phylo
td <- as_tibble(mcc_tree)
n_tips <- Ntip(tree)

## Step 1: pull internal nodes above a posterior threshold
post_thresh <- 0.9
target_k <- 8  # target number of clades (adjust to taste)

candidates <- td %>% filter(node > n_tips, posterior >= post_thresh)

## Step 2: get descendant tip sets for each candidate node
desc_tips <- Descendants(tree, candidates$node, type = "tips")
candidates$n_tips <- sapply(desc_tips, length)
candidates$tips   <- desc_tips

## Step 3: filter to a reasonable clade size range
min_size <- floor(n_tips / (target_k * 3))   # ~149 tips
max_size <- ceiling(n_tips / 2)
candidates <- candidates %>% filter(n_tips >= min_size, n_tips <= max_size) %>%
  arrange(desc(posterior), desc(n_tips))      # prefer high support, then large

## Step 4: greedy non-overlapping selection
selected <- list()
covered  <- integer(0)

for (i in seq_len(nrow(candidates))) {
  tips_i <- candidates$tips[[i]]
  if (length(intersect(tips_i, covered)) == 0) {
    selected[[length(selected) + 1]] <- list(
      node      = candidates$node[i],
      posterior = candidates$posterior[i],
      n_tips    = candidates$n_tips[i],
      tips      = tips_i)
    covered <- union(covered, tips_i)
  }
  # stop once we've reached the target count and good coverage
  if (length(selected) >= target_k & length(covered) >= n_tips * 0.90) break
}

## Summary of what was selected and level of coverage
cat("Clades selected:", length(selected), "\n")
cat("Tips covered:", length(covered), "of", n_tips,
    paste0("(", round(100 * length(covered) / n_tips, 1), "%)\n"))
for (j in seq_along(selected)) {
  cat(sprintf("  Clade %d: node %d, n=%d, posterior=%.3f\n",
              j, selected[[j]]$node, selected[[j]]$n_tips, selected[[j]]$posterior))
}


## Step 5: if coverage is low, relax the threshold and fill gaps
uncovered <- setdiff(seq_len(n_tips), covered)

if (length(uncovered) > n_tips * 0.05) {
  cat("\nRelaxing to posterior >= 0.25 for remaining??", length(uncovered), "tips...\n")
  
  cands2 <- td %>%
    filter(node > n_tips, posterior >= 0.5, !(node %in% sapply(selected, `[[`, "node")))
  
  desc2 <- Descendants(tree, cands2$node, type = "tips")
  cands2$n_tips <- sapply(desc2, length)
  cands2$tips   <- desc2
  
  # only keep nodes whose tips are entirely within uncovered set
  cands2 <- cands2 %>% filter(sapply(tips, function(t) all(t %in% uncovered)),
           n_tips >= min_size) %>% arrange(desc(n_tips))
  
  for (i in seq_len(nrow(cands2))) {
    tips_i <- cands2$tips[[i]]
    if (length(intersect(tips_i, covered)) == 0) {
      selected[[length(selected) + 1]] <- list(
        node      = cands2$node[i],
        posterior = cands2$posterior[i],
        n_tips    = cands2$n_tips[i],
        tips      = tips_i)
      covered <- union(covered, tips_i)
    }
    if (length(setdiff(seq_len(n_tips), covered)) < 10) break
  }
  cat("After relaxation:", length(selected), "clades,",
      length(covered), "tips covered\n")
}

## Step 6: build clade assignment vector and merge into tip_meta
tip_labels <- tree$tip.label
clade_vec <- rep(NA_character_, n_tips)

for (j in seq_along(selected)) {
  clade_vec[selected[[j]]$tips] <- paste0("Clade ", j)
}

clade_df <- data.frame(label = tip_labels, clade = clade_vec, stringsAsFactors = FALSE)

# After building clade_df, renumber by median y-position
plot_data <- ggtree(mcc_tree)$data %>% filter(isTip)
clade_df <- clade_df %>% left_join(plot_data %>% select(label, y), by = "label") %>%
  group_by(clade) %>% mutate(median_y = median(y, na.rm = TRUE)) %>%
  ungroup() %>% mutate(clade = paste0("Clade ", dense_rank(median_y))) %>%
  select(label, clade)


## Step 7: assign orphans to the clade of their nearest tip on the tree
uncovered <- setdiff(seq_len(n_tips), covered)

if (length(uncovered) > 0) {
  cat("Assigning", length(uncovered), "orphan tips to nearest clade...\n")
  
  # Build a full distance matrix (can take a moment with ~2700 tips)
  all_dist <- dist.nodes(tree)
  
  # For each orphan, find the nearest tip that IS in a clade
  orphan_clade <- sapply(uncovered, function(tip) {
    # distances from this orphan to all covered tips
    dists_to_covered <- all_dist[tip, setdiff(seq_len(n_tips), uncovered)]
    nearest_tip <- setdiff(seq_len(n_tips), uncovered)[which.min(dists_to_covered)]
    # return whichever clade that nearest tip belongs to
    clade_vec[nearest_tip]
  })
  
  # Assign
  clade_vec[uncovered] <- orphan_clade
}

clade_df <- data.frame(label = tip_labels, clade = clade_vec, stringsAsFactors = FALSE)
tip_meta <- left_join(tip_meta, clade_df, by="label")


## Step 8: plot
baseline_clades <- ggtree(mcc_tree) %<+% tip_meta +
  geom_tippoint(aes(color = clade), size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set3") +
  theme_tree2() + coord_cartesian(xlim = c(0, 26)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm"))
baseline_clades


## Step 9: export node and clade summary
# Count how many tips were directly assigned vs orphan-assigned per clade
direct_df <- data.frame(tip_index = unlist(lapply(selected, `[[`, "tips")),
  stringsAsFactors = FALSE)
direct_set <- direct_df$tip_index

clade_summary <- data.frame(
  clade     = paste0("Clade ", seq_along(selected)),
  node      = sapply(selected, `[[`, "node"),
  posterior = sapply(selected, `[[`, "posterior"),
  n_direct  = sapply(selected, `[[`, "n_tips"),
  stringsAsFactors = FALSE)

# Renumber clades to match the y-position renumbering from Step 6
clade_summary <- clade_summary %>%
  left_join(
    clade_df %>% filter(label %in% tip_labels[unlist(lapply(selected, `[[`, "tips"))]) %>%
      distinct(clade), by = character())

# Simpler approach: count from the final clade_df directly
clade_summary <- clade_df %>%
  mutate(is_orphan = !(match(label, tip_labels) %in% direct_set)) %>%
  group_by(clade) %>% summarise(
    n_total   = n(),
    n_direct  = sum(!is_orphan),
    n_orphan  = sum(is_orphan),
    .groups   = "drop") %>% arrange(clade)

# Add node and posterior info (map original clade numbers to renumbered ones)
node_info <- data.frame(
  original  = paste0("Clade ", seq_along(selected)),
  node      = sapply(selected, `[[`, "node"),
  posterior = round(sapply(selected, `[[`, "posterior"), 3),
  stringsAsFactors = FALSE)

# Match by direct tip count to link original selected clades to renumbered clades
node_lookup <- data.frame(
  label = tip_labels[unlist(lapply(seq_along(selected), function(j) selected[[j]]$tips[1]))],
  original = paste0("Clade ", seq_along(selected)),
  stringsAsFactors = FALSE) %>% left_join(clade_df, by = "label") %>%
  select(original, clade) %>% left_join(node_info, by = "original")

clade_summary <- clade_summary %>%
  left_join(node_lookup %>% select(clade, node, posterior), by = "clade")

print(clade_summary, row.names = FALSE)
cat("\nTotal tips:", sum(clade_summary$n_total), "of", n_tips, "\n")

