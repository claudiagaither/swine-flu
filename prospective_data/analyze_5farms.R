# =============================================================================
# MultiTypeTree Migration Analysis - 5 FARMS (Backward Only)
# 1FS and FS1 merged into "McCook"
# =============================================================================

library(coda)
library(ggplot2)
library(reshape2)
library(igraph)
library(viridis)
library(dplyr)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

log_file <- "C:/Users/cgait/OneDrive/Desktop/swine flu/multitypetree/h3_MTT_5farms_backward_25M.log"
trace <- read.table(log_file, header = TRUE)

# Remove burn-in (first 10% of samples)
burnin <- round(0.25 * nrow(trace))
trace_post_burnin <- trace[(burnin + 1):nrow(trace), ]

cat("=== CONVERGENCE DIAGNOSTICS ===\n\n")

# Check ESS for key parameters
ess_values <- effectiveSize(trace_post_burnin[, 2:ncol(trace_post_burnin)])

# Migration rate ESS
mig_ess <- ess_values[grep("rateMatrix", names(ess_values))]
cat("Migration rate ESS summary:\n")
cat("  Min:", min(mig_ess), "\n")
cat("  Median:", median(mig_ess), "\n")
cat("  Max:", max(mig_ess), "\n")
cat("  Number with ESS < 200:", sum(mig_ess < 200), "out of", length(mig_ess), "\n\n")

# Population size ESS
pop_ess <- ess_values[grep("popSize_", names(ess_values))]
cat("Population size ESS:\n")
print(pop_ess)
cat("\n")

# =============================================================================
# 2. FARM METADATA (5 FARMS)
# =============================================================================

farm_metadata <- data.frame(
  Farm_ID = c("McCook", "AH1", "GD1", "KH1", "RF1"),
  County = c("McCook", "Pipestone", "Scott", "Osceola", "Turner"),
  State = c("SD", "MN", "MO", "IA", "SD"),
  N_Samples = c(12, 5, 7, 1, 12),
  stringsAsFactors = FALSE
)

farm_metadata$Location <- paste0(farm_metadata$County, ", ", farm_metadata$State)

cat("=== FARM METADATA ===\n\n")
print(farm_metadata)
cat("\n")

# =============================================================================
# 3. EXTRACT BACKWARD MIGRATION RATES
# =============================================================================

cat("=== MIGRATION RATE ANALYSIS ===\n\n")

farm_ids <- c("McCook", "AH1", "GD1", "KH1", "RF1")
n_farms <- length(farm_ids)

# Only backward rates (no forward rates in this model)
backward_columns <- grep("rateMatrix_backward", colnames(trace_post_burnin), value = TRUE)

cat("Found", length(backward_columns), "backward rate parameters\n")
cat("Expected:", n_farms * (n_farms - 1), "rates for", n_farms, "farms\n\n")

# Extract backward rates
backward_rates <- apply(trace_post_burnin[, backward_columns], 2, mean)

# Parse parameter names
parse_migration_param <- function(param_name) {
  parts <- strsplit(param_name, "_")[[1]]
  n_parts <- length(parts)
  from_farm <- parts[n_parts - 1]
  to_farm <- parts[n_parts]
  return(c(from = from_farm, to = to_farm))
}

# Create migration matrix
migration_matrix <- matrix(0, nrow = n_farms, ncol = n_farms,
                          dimnames = list(farm_ids, farm_ids))

for (i in 1:length(backward_rates)) {
  param_name <- names(backward_rates)[i]
  farms <- parse_migration_param(param_name)
  from_idx <- which(farm_ids == farms["from"])
  to_idx <- which(farm_ids == farms["to"])
  
  if (length(from_idx) > 0 && length(to_idx) > 0) {
    migration_matrix[from_idx, to_idx] <- backward_rates[i]
  }
}

cat("BACKWARD MIGRATION RATE MATRIX (events per unit time):\n")
cat("(Rows = source, Columns = destination)\n\n")
print(round(migration_matrix, 3))
cat("\n")

# =============================================================================
# 4. MIGRATION SUMMARY
# =============================================================================

cat("=== MIGRATION SUMMARY ===\n\n")

migration_summary <- data.frame(
  Farm_ID = farm_ids,
  Total_Emigrations = rowSums(migration_matrix),
  Total_Immigrations = colSums(migration_matrix),
  Net_Migration = colSums(migration_matrix) - rowSums(migration_matrix),
  stringsAsFactors = FALSE
)

# Add metadata by matching Farm_ID
migration_summary$Location <- farm_metadata$Location[match(migration_summary$Farm_ID, farm_metadata$Farm_ID)]
migration_summary$N_Samples <- farm_metadata$N_Samples[match(migration_summary$Farm_ID, farm_metadata$Farm_ID)]

# Reorder columns
migration_summary <- migration_summary[, c("Farm_ID", "Location", "N_Samples", 
                                          "Total_Emigrations", "Total_Immigrations", "Net_Migration")]

cat("Migration summary by farm:\n")
print(migration_summary)
cat("\n")

# =============================================================================
# 5. VISUALIZATIONS
# =============================================================================

cat("=== CREATING VISUALIZATIONS ===\n\n")

# Prepare data for plotting
migration_df <- melt(migration_matrix)
colnames(migration_df) <- c("From", "To", "Rate")
migration_df$From <- as.character(migration_df$From)
migration_df$To <- as.character(migration_df$To)
migration_df <- migration_df[migration_df$From != migration_df$To, ]

# Add metadata
migration_df <- merge(migration_df, farm_metadata[, c("Farm_ID", "County", "Location")], 
                     by.x = "From", by.y = "Farm_ID", all.x = TRUE)
colnames(migration_df)[4:5] <- c("From_County", "From_Location")

migration_df <- merge(migration_df, farm_metadata[, c("Farm_ID", "County", "Location")], 
                     by.x = "To", by.y = "Farm_ID", all.x = TRUE)
colnames(migration_df)[6:7] <- c("To_County", "To_Location")

# Check for same-county pairs (none expected now since each farm is in different county)
migration_df$Same_County <- migration_df$From_County == migration_df$To_County
migration_df$Same_State <- substr(migration_df$From_Location, 
                                  nchar(migration_df$From_Location)-1, 
                                  nchar(migration_df$From_Location)) == 
                           substr(migration_df$To_Location, 
                                  nchar(migration_df$To_Location)-1, 
                                  nchar(migration_df$To_Location))

migration_df <- migration_df %>% filter(From_Location != "RF1")

# -----------------------------------------------------------------------------
# PLOT 1: Migration Rate Heatmap
# -----------------------------------------------------------------------------

p1 <- ggplot(migration_df, aes(x = To_Location, y = From_Location, fill = Rate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(Rate, 2)), color = "white", size = 4, fontface = "bold") +
  # Highlight same-state pairs
  geom_tile(data = migration_df[migration_df$Same_State & !migration_df$Same_County, ], 
            aes(x = To, y = From), 
            color = "yellow", fill = NA, size = 1) +
  scale_fill_viridis(option = "plasma", name = "Migration\nRate") +
  labs(title = "Migration Rates Between 5 Farms",
       subtitle = "Yellow border = same state, different county",
       x = "Destination Farm",
       y = "Source Farm") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

print(p1)
ggsave("migration_heatmap_5farms.png", p1, width = 8, height = 7, dpi = 300)
cat("Saved: migration_heatmap_5farms.png\n")

# -----------------------------------------------------------------------------
# PLOT 2: Network Diagram
# -----------------------------------------------------------------------------

# Filter for rates > quantile for cleaner visualization
threshold <- quantile(migration_df$Rate, 0.5)
edges_filtered <- migration_df[migration_df$Rate > threshold, ]

# Create igraph
g <- graph_from_data_frame(edges_filtered[, c("From", "To", "Rate")], 
                           directed = TRUE, 
                           vertices = farm_metadata)

# Set vertex attributes
V(g)$county <- farm_metadata$County[match(V(g)$name, farm_metadata$Farm_ID)]
V(g)$location <- farm_metadata$Location[match(V(g)$name, farm_metadata$Farm_ID)]

# Color by county
county_colors <- c("McCook" = "#E41A1C",   # Red
                   "Pipestone" = "#377EB8", # Blue
                   "Scott" = "#4DAF4A",     # Green
                   "Osceola" = "#984EA3",   # Purple
                   "Turner" = "#FF7F00")    # Orange

V(g)$color <- county_colors[V(g)$county]

# Edge attributes
E(g)$width <- scales::rescale(E(g)$Rate, to = c(0.5, 5))
E(g)$color <- "gray50"

# Plot
png("migration_network_5farms.png", width = 10, height = 10, units = "in", res = 300)
par(mar = c(1, 1, 3, 1))
plot(g, 
     vertex.size = 40,
     vertex.label.color = "white",
     vertex.label.font = 2,
     vertex.label.cex = 1.2,
     edge.arrow.size = 0.5,
     edge.curved = 0.3,
     layout = layout_with_fr(g),
     main = "Migration Network (5 Farms)\nEdge width ∝ migration rate")

legend("bottomright", 
       legend = names(county_colors), 
       fill = county_colors, 
       title = "County",
       bty = "n",
       cex = 1.2)
dev.off()
cat("Saved: migration_network_5farms.png\n")

# -----------------------------------------------------------------------------
# PLOT 3: Top Migration Routes
# -----------------------------------------------------------------------------

top_routes <- migration_df[order(migration_df$Rate, decreasing = TRUE), ][1:10, ]
top_routes$Route <- paste(top_routes$From, "→", top_routes$To)
top_routes$Route <- factor(top_routes$Route, levels = rev(top_routes$Route))

p3 <- ggplot(top_routes, aes(x = Rate, y = Route, fill = Same_State)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral"),
                    labels = c("Different state", "Same state"),
                    name = "") +
  labs(title = "Top 10 Migration Routes",
       x = "Migration Rate (events per unit time)",
       y = "Migration Route") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom")

print(p3)
ggsave("top_routes_5farms.png", p3, width = 10, height = 6, dpi = 300)
cat("Saved: top_routes_5farms.png\n")

# -----------------------------------------------------------------------------
# PLOT 4: State-Level Summary
# -----------------------------------------------------------------------------

# Extract state from location
migration_df$From_State <- sapply(strsplit(migration_df$From_Location, ", "), function(x) x[2])
migration_df$To_State <- sapply(strsplit(migration_df$To_Location, ", "), function(x) x[2])

state_migration <- migration_df %>%
  group_by(From_State, To_State) %>%
  summarise(
    Total_Rate = sum(Rate),
    N_Routes = n(),
    Avg_Rate = mean(Rate),
    .groups = "drop"
  ) %>%
  filter(From_State != To_State)

p4 <- ggplot(state_migration, aes(x = To_State, y = From_State, fill = Avg_Rate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f\n(%d)", Avg_Rate, N_Routes)), 
            color = "white", size = 4, fontface = "bold") +
  scale_fill_viridis(option = "magma", name = "Average\nRate") +
  labs(title = "Average Migration Rates Between States",
       subtitle = "Numbers show: Average rate (number of farm pairs)",
       x = "Destination State",
       y = "Source State") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p4)
ggsave("state_level_migration.png", p4, width = 8, height = 6, dpi = 300)
cat("Saved: state_level_migration.png\n\n")

# =============================================================================
# 6. EXPORT RESULTS
# =============================================================================

write.csv(migration_matrix, "migration_matrix_5farms.csv", row.names = TRUE)
write.csv(migration_summary, "migration_summary_5farms.csv", row.names = FALSE)
write.csv(migration_df, "pairwise_migration_5farms.csv", row.names = FALSE)

cat("=== RESULTS EXPORTED ===\n\n")
cat("CSV files:\n")
cat("  1. migration_matrix_5farms.csv\n")
cat("  2. migration_summary_5farms.csv\n")
cat("  3. pairwise_migration_5farms.csv\n\n")

cat("PNG files:\n")
cat("  1. migration_heatmap_5farms.png\n")
cat("  2. migration_network_5farms.png\n")
cat("  3. top_routes_5farms.png\n")
cat("  4. state_level_migration.png\n\n")

# =============================================================================
# 7. SUMMARY STATISTICS
# =============================================================================

cat("=== FINAL SUMMARY ===\n\n")

cat("Farm-level statistics:\n")
cat("  Highest emigration:", migration_summary$Farm_ID[which.max(migration_summary$Total_Emigrations)], 
    "(", round(max(migration_summary$Total_Emigrations), 2), "events/time)\n")
cat("  Highest immigration:", migration_summary$Farm_ID[which.max(migration_summary$Total_Immigrations)], 
    "(", round(max(migration_summary$Total_Immigrations), 2), "events/time)\n")
cat("  Strongest source (net emigration):", 
    migration_summary$Farm_ID[which.min(migration_summary$Net_Migration)],
    "(net:", round(min(migration_summary$Net_Migration), 2), ")\n")
cat("  Strongest sink (net immigration):", 
    migration_summary$Farm_ID[which.max(migration_summary$Net_Migration)],
    "(net:", round(max(migration_summary$Net_Migration), 2), ")\n\n")

cat("Top migration route:", top_routes$Route[1], 
    "with rate", round(top_routes$Rate[1], 2), "\n\n")

cat("State-level patterns:\n")
if (nrow(state_migration) > 0) {
  top_state <- state_migration[which.max(state_migration$Avg_Rate), ]
  cat("  Highest state-to-state migration:", 
      top_state$From_State, "→", top_state$To_State, 
      "(avg rate:", round(top_state$Avg_Rate, 2), ")\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
