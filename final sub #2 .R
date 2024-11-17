# Load required packages
library(tidyverse)
library(rentrez)
library(ape)
library(msa)
library(phangorn)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(knitr)
library(geosphere)
library(raster)
library(sp)
library(readr) # For read_tsv()
library(seqinr)
library(Biostrings)
library(ggtree)
library(DECIPHER)

# Data Loading and Initial Filtering -------------------------------------------
#_ BOLD Data Filtering ----

# Initializing standard requirement variables for sequence data during filtering
length.var <- 50
missing.data <- 0.01

# Load and filter data from BOLD dataset, and create a subset of unique species.
# Original call to pull BOLD data was not provided.
dflizard <- read_tsv("./bold_data-2.tsv", show_col_types = FALSE) %>%
  as.data.frame()

# This step filters sequence length and variance to establish a strong phylogenetic tree in downstream analysis. Specific data are also filtered for (species name, location data and sequence) and filtered out (any entry that is missing data).
dfBOLD <- dflizard %>%
  filter(markercode == "COI-5P") %>%
  filter(!is.na(species_name), !is.na(processid), !is.na(lat), !is.na(lon), !is.na(nucleotides), !is.na(country)) %>%
  mutate(COI_Sequence = str_remove_all(nucleotides, "^N+|N+$|-")) %>%
  filter(str_count(COI_Sequence, "N") <= (missing.data * str_count(nucleotides))) %>%
  filter(str_count(COI_Sequence) >= median(str_count(COI_Sequence)) - length.var & str_count(COI_Sequence) <= median(str_count(COI_Sequence)) + length.var) %>%
  mutate(header = paste(processid, species_name, sep = "_")) %>%
  dplyr::select(header, processid, species_name, country, lat, lon, nucleotides, markercode) %>% 
  group_by(species_name) %>%
  sample_n(1)

# Check filtered data
view(dfBOLD)

#_ Fetch COI Gene Sequences from NCBI -------------------------------------------
# All of the below is commented out from the original script.
# Fetch sequence data from NCBI for Phelsuma species to build a phylogenetic tree. Saves sequences in FASTA format for reproducibility and future use.
# Search for Phelsuma COI gene sequences
# search <- entrez_search(db = "nucleotide", term = "Phelsuma COI", retmax = 100)
# 
# # Proceed if there are results
# if (search$count > 0) {
#   
#   # Retrieve metadata for the fetched sequence IDs
#   metadata <- entrez_summary(db = "nucleotide", id = search$ids)
#   
#   # Extract organism names from the metadata
#   species_info <- sapply(metadata, function(x) x$organism)  # Use the 'organism' field for species names
#   
#   # Create a data frame with species names
#   sequence_data <- data.frame(Species = species_info, stringsAsFactors = FALSE)
#   
#   # View the combined data
#   print(sequence_data)
#   
#   # saving them as Fasta files 
#   
#   # Retrieve the COI sequences in FASTA format using sequence IDs
#   fasta_sequences <- entrez_fetch(db = "nucleotide", id = search$ids, rettype = "fasta", retmode = "text")
#   
#   # Save the sequences to a FASTA file on your computer
#   write(fasta_sequences, file = "../6210 assignment 2/phelsuma_sequences.fasta")
#   
#   # View the first few lines of the saved FASTA file (optional)
#   cat(readLines("../6210 assignment 2/phelsuma_sequences.fasta", n = 10), sep = "\n")
#   
#   cat("FASTA sequences saved successfully!")
# }

#_ NCBI Data Filtering ----
# Import FASTA data from file, and create a subset of unique species.
NCBI.string.set <- readDNAStringSet("./phelsuma_sequences.fasta")
dfNCBI <- data.frame(COI_Title = names(NCBI.string.set), COI_Seq = paste(NCBI.string.set))
dfNCBI$Species_Name <- word(dfNCBI$COI_Title, start = 2L, end = 3L)

# This step filters sequence length and variance to establish a strong phylogenetic tree in downstream analysis.
dfNCBI <- dfNCBI %>%
  mutate(COI_Sequence = str_remove_all(COI_Seq, "^N+|N+$|-")) %>%
  filter(str_count(COI_Sequence, "N") <= (missing.data * str_count(COI_Seq))) %>%
  filter(str_count(COI_Sequence) >= median(str_count(COI_Sequence)) - length.var & str_count(COI_Sequence) <= median(str_count(COI_Sequence)) + length.var) %>%
  filter(!str_detect(Species_Name, "UNVERIFIED")) %>%
  dplyr::select(Species_Name, COI_Sequence) %>% 
  group_by(Species_Name) %>%
  sample_n(1)

#_ Combining Data from BOLD & NCBI to streamline future analysis.
# Combine data as a dataframe
combined_data <- dfNCBI %>% 
  left_join(dfBOLD, by = c("Species_Name" = "species_name"))

# Write combined data to FASTA format
write.fasta(sequences = as.list(combined_data$nucleotides), names = combined_data$Species_Name, file.out = "combined.fasta")

# Phylogenetic Analysis --------------------------------------------------------
# Keep the original names of the fasta file, as the MSA method removes them.
original_names <- names(readDNAStringSet("combined.fasta"))

# Perform the Multiple Sequence Alignment, and convert format for future analysis
alignment <- msa(readDNAStringSet("combined.fasta"), type = "DNA")
aligned_sequences <- as(alignment, "DNAStringSet")
names(aligned_sequences) <- original_names # Reapply the names
print(names(aligned_sequences)) #Verify the names.

# View the alignnment to visually check for errors, then write it to file.
BrowseSeqs(aligned_sequences)
writeXStringSet(aligned_sequences, filepath = "aligned_combined.fasta")

# Read the aligned sequences to the correct format for making a phylogenetic tree
alignment_data <- read.phyDat("aligned_combined.fasta", format = "fasta")

# Create a distance matrix using Maximum Likelihood, and construct an initial tree using Neighbor Joining
dist_matrix <- dist.ml(alignment_data)
nj_tree <- nj(dist_matrix)

# Fit the initial tree using the ML approach
fit <- pml(nj_tree, data = alignment_data)

# Test the evolutionary model used to find the best match for the data
model_test_results <- modelTest(fit, model = "all")

# Print the best model based on AIC score.
best_model <- model_test_results[which.min(model_test_results$AIC), ]
print(best_model)

# Fit the best model to the data and convert the data to the correct format for visualization.
fit_model <- update(fit, model = "GTR", k = 4, inv = 0.2)
tree_phylo <- as.phylo(fit_model$tree)

# Plot the tree using a simple layout.
p <- ggtree(tree_phylo) +
  ggtitle("Maximum Likelihood Tree of Phelsuma Species") +
  theme_tree() +
  geom_treescale(x = 0.2, y = 0.5, label = "Branch Length") +
  geom_tiplab(size = 2)
print(p)

# This tree shows a visual indication of sister species. Now we will identify them using the Most Recent Common Ancestor (MRCA) method.

# Sister Species Identification ------------------------------------------------
# Generate all potential pairs of species from the ML tree labels.
species_pairs <- combn(tree_phylo$tip.label, 2, simplify = FALSE)

# Initialize lists that will store the sister species and their MRCA
sister_species_list <- list()
mrca_nodes <- list()

# Use a for loop to identify sister species based on the indices of the pair in the tree's tip labels.
for (pair in species_pairs) {
  tip_indices <- match(pair, tree_phylo$tip.label)
  
  # Store the MRCA for the identified pair.
  mrca <- getMRCA(tree_phylo, tip_indices)
  
  # Find the descendants of the node.
  descendants <- tree_phylo$edge[tree_phylo$edge[, 1] == mrca, 2]

  # This is if statement ensures that there is a pair of sister species, and saves their names in the sister_species list, along with their MRCA.
  if (length(descendants) == 2) {
    sister_species <- tree_phylo$tip.label[descendants]
    sister_species_list[[paste(pair, collapse = " & ")]] <- sister_species
    mrca_nodes[[paste(pair, collapse = " & ")]] <- mrca
  } else {
    next    # Skips if there are not a pair of descendents.
  }
}

# Stores the sister_species as pairs
sister_species_pairs <- data.frame(
  species_1 = unlist(lapply(sister_species_list, function(x) x[1])),
  species_2 = unlist(lapply(sister_species_list, function(x) x[2]))
)

# Filters the sisters based on geographic data.
geographic_sisters <- combined_data %>%
  filter(Species_Name %in% c(sister_species_pairs$species_1, sister_species_pairs$species_2)) %>%
  dplyr::select(Species_Name, country, lat, lon) %>%
  filter(!is.na(lon) & !is.na(lat))

# This adds a column for pair IDs, which will be used to help detect the presence of pairs and their geographical relationship later on.
pair_ids <- rep(1:nrow(sister_species_pairs), each = 2)
geographic_sisters$pair_id <- NA

# This for loop assigns pair IDs based on the filtered geographic data for matching sister pairs.
for (i in 1:nrow(sister_species_pairs)) {
  species_1 <- sister_species_pairs$species_1[i]
  species_2 <- sister_species_pairs$species_2[i]
  
  # Finds the index of each species
  species_1_index <- which(geographic_sisters$Species_Name == species_1)
  species_2_index <- which(geographic_sisters$Species_Name == species_2)

  # Assigns the pair ID to matching sisters.
  geographic_sisters$pair_id[species_1_index] <- i
  geographic_sisters$pair_id[species_2_index] <- i
}

# Check that the Pair IDs were added.
view(geographic_sisters)

# Categorize Sister Species by Region Overlap-----------------------------------
# This adds the location information to each pair, and also highlights which of the species do not have a sister species.
geographic_sisters <- geographic_sisters %>%
  group_by(pair_id) %>%
  mutate(location_status = case_when(
    n() == 2 & n_distinct(country) == 1 ~ "Same Location",
    n() == 2 & n_distinct(country) > 1 ~ "Different Location",
    n() == 1 ~ "No Pair",
    TRUE ~ "Other"
  )) %>%
  ungroup()

# Create maps for different location statuses ----
# For each location status, I create a map to examine the geographical spread. Three seperate maps are created to avoid crowding of labels.

# This function creates a base world map for use in each of the three maps.
world_map <- map_data("world")

create_location_map <- function(data, color, title, world_map) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
    geom_point(data = data, aes(x = lon, y = lat), color = color, size = 3) +
    geom_text_repel(
      data = data, aes(x = lon, y = lat, label = Species_Name),
      size = 3, box.padding = 0.3, point.padding = 0.2, segment.color = "grey",
      segment.size = 0.5, max.overlaps = 50, direction = "both"
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# This creates a map of all species that have a sister in the same geographical area.
same_location_map <- create_location_map(
  data = geographic_sisters %>% filter(location_status == "Same Location"),
  color = "green",
  title = "Sister Species in the Same Location",
  world_map = world_map
)

# This creates a map of all species that have a sister species in a different geographical area.
different_location_map <- create_location_map(
  data = geographic_sisters %>% filter(location_status == "Different Location"),
  color = "red",
  title = "Sister Species in Different Locations",
  world_map = world_map
)

# This map highlights only species that do not have a sister.
no_pair_map <- create_location_map(
  data = geographic_sisters %>% filter(location_status == "No Pair"),
  color = "gray",
  title = "Species Without Pairs",
  world_map = world_map
)

# Print the maps
print(same_location_map)
print(different_location_map)
print(no_pair_map)


# Interpretation: These plots helps reveal whether sister species occupy the same or different regions, providing insights into potential allopatric or sympatric speciation.

# Creaing a Table to Represent Location Status ----
# Filter to show only pairs with complete region data
sister_summary_table_complete <- geographic_sisters %>%
  filter(location_status != "No Pair") %>%
  dplyr::select(Species_Name, country, pair_id, location_status) %>%
  group_by(pair_id)

# Print the filtered table
print(sister_summary_table_complete)

# Display the table in a clean format
kable(sister_summary_table_complete,
  col.names = c("Species Name", "Country", "Pair ID", "Location Status"),
  caption = "Geographic Overlap of Sister Species in Phelsuma",
  align = "l") # left-align columns

# Quantify Region Overlap--------------------------------------------------------------
# Summarize and categorize sister species based on geographic overlap to quantify speciation patterns.
# Summary of region overlap counts
overlap_summary <- sister_summary_table_complete %>%
  group_by(location_status) %>%
  summarize(count = n())

# Print the summary
print(overlap_summary)

# make it better looking and adding a visualization
overlap_summary <- as.data.frame(overlap_summary)
colnames(overlap_summary) <- c("Geographic Overlap", "Number of Pairs")
print(overlap_summary)

# Interpretation: The table quantifies the occurrence of speciation in shared vs. different regions, addressing the main question of geographic influence on speciation.








## calculate geographical distances-------------------------------------------------------
# TODO Remaining script is untouched, and therefore will throw errors. Updates to the methods below must be completed with new dataframe structures taken into account.

# Calculate distances between sister species and categorize them to test spatial influence on speciation.
sister_lines_data_clean <- sister_lines_data_clean %>%
  rowwise() %>%
  mutate(distance_km = distHaversine(
    c(lon_1, lat_1),
    c(lon_2, lat_2)
  )) %>%
  ungroup()

# View the updated data with distances
print(sister_lines_data_clean)

# Interpret Specific Sister Species Pairs------------------------------------------------------
# Classify distances as "close" or "far" based on a 100 km threshold
sister_lines_data_clean <- sister_lines_data_clean %>%
  mutate(distance_category = if_else(distance_km <= 100, "close", "far"))

# View the updated dataset with the new distance category
print(sister_lines_data_clean)

# Summarize the count of "close" and "far" pairs
distance_summary <- sister_lines_data_clean %>%
  group_by(distance_category) %>%
  summarize(count = n())

print(distance_summary)

# visualize
ggplot(sister_lines_data_clean, aes(x = distance_category, fill = country_1 == country_2)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Geographic Distance Classification of Sister Species in Phelsuma",
    x = "Distance Category",
    y = "Number of Pairs",
    fill = "Same Region"
  ) +
  scale_fill_manual(values = c("TRUE" = "salmon", "FALSE" = "steelblue")) +
  theme_minimal()

# Plot the geographic distribution and connections between sister species pairs with improved readability

# Adjust plot dimensions for improved readability
options(repr.plot.width = 14, repr.plot.height = 9)

# Generate the improved plot
ggplot(data = sister_lines_data_clean, aes(x = lon_1, y = lat_1)) +
  borders("world", colour = "grey85", fill = "grey95") +
  geom_segment(aes(
    xend = lon_2, yend = lat_2, color = distance_category,
    linetype = distance_category, size = distance_category
  ), alpha = 0.5) + # Reduced alpha for lines
  scale_color_manual(values = c("close" = "darkred", "far" = "blue")) +
  scale_linetype_manual(values = c("close" = "solid", "far" = "dotted")) +
  scale_size_manual(values = c("close" = 1.2, "far" = 0.8)) +
  geom_point(aes(color = country_1), size = 3, alpha = 0.7) +
  geom_point(aes(x = lon_2, y = lat_2, color = country_2), size = 3, alpha = 0.7) +
  geom_text_repel(
    aes(label = species_1),
    size = 2, fontface = "italic", max.overlaps = Inf,
    box.padding = 1, segment.color = "grey50", force = 2, # Increased box.padding for blue dot area
    bg.color = "white", bg.r = 0.15
  ) +
  geom_text_repel(
    aes(x = lon_2, y = lat_2, label = species_2),
    size = 2, fontface = "italic", max.overlaps = Inf,
    box.padding = 1, segment.color = "grey50", force = 2,
    bg.color = "white", bg.r = 0.15
  ) +
  # Adding repel feature for distance labels near blue dot
  geom_text_repel(
    aes(x = (lon_1 + lon_2) / 2, y = (lat_1 + lat_2) / 2, label = paste0(round(distance_km, 1), " km")),
    color = "darkgreen", size = 2.2, fontface = "bold", hjust = 0.5, alpha = 0.75,
    max.overlaps = Inf, segment.size = 0.3, box.padding = 0.4, force = 1 # Adjusted for readability
  ) +
  labs(
    title = "Geographic Distances Between Sister Species in Phelsuma",
    x = "Longitude", y = "Latitude", color = "Distance Category"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("Comoros" = "salmon", "France" = "green", "Madagascar" = "skyblue", "United States" = "purple"))

# to bring everything together and give definitive conclusion-------------------------------------
# quantify overall findings
# Count total pairs and pairs in each region overlap category
total_pairs <- nrow(sister_summary_table_complete)
same_region_pairs <- nrow(sister_summary_table_complete %>% filter(region_overlap == "Same Region"))
different_region_pairs <- nrow(sister_summary_table_complete %>% filter(region_overlap == "Different Region"))

# Calculate percentages
same_region_percentage <- (same_region_pairs / total_pairs) * 100
different_region_percentage <- (different_region_pairs / total_pairs) * 100

# Display the results
cat("Same Region Pairs:", same_region_percentage, "%\n")
cat("Different Region Pairs:", different_region_percentage, "%\n")

# check patterns in my data----------------------------------------------------------
# 1. Test if Sister Species are More Likely in the Same or Different Regions
# Count table for "Same Region" and "Different Region"
region_table <- table(sister_lines_data_clean$distance_category)

# Run the Chi-square test
chi_test <- chisq.test(region_table)

# View the test results
print(chi_test)

# Interpretation: The Chi-square test evaluates whether sister species pairs are more likely
# to exist in the same or different regions, addressing the main hypothesis on speciation.The p-value (0.2) was insignificant suggesting the observed distribution of sister species pairs across regions could be due to chance rather than a clear preference for same or different regions

# adding and exploring using the climate data-----------------------------------------------------------------

# Define the path where the .tif file is located after extraction
bio1_path <- "../6210 assignment 2/wc2.1_2.5m_bio_1/wc2.1_2.5m_bio_1.tif"

# Load the raster
bio1 <- raster(bio1_path)

# precipitation info
# Define the path where the .tif file is located after extraction
bio12_path <- "../6210 assignment 2/wc2.1_2.5m_bio_12/wc2.1_2.5m_bio_12.tif"

# Load the raster
bio12 <- raster(bio12_path)

# climate analysis-------------------------------------------------------------------------
# Code to Create sister_species_coords
# Get a list of unique sister species names from `sister_pairs_df`
sister_species <- unique(c(sister_pairs_df$species_1, sister_pairs_df$species_2))

#  Filter `dflizardfilter` to include only sister species
sister_species_coords <- dflizardfilter %>%
  filter(species_name %in% sister_species) %>%
  dplyr::select(species_name, lat, lon) %>%
  dplyr::rename(latitude = lat, longitude = lon)

# View the `sister_species_coords` data frame to confirm
print(sister_species_coords)

# visualization
# Extract climate data for each location in `sister_species_coords`
temperature_values <- raster::extract(bio1, sister_species_coords[, c("longitude", "latitude")])
precipitation_values <- raster::extract(bio12, sister_species_coords[, c("longitude", "latitude")])

# Add extracted values to `sister_species_coords`
sister_species_coords$temperature <- temperature_values
sister_species_coords$precipitation <- precipitation_values

# Quick Visualization of Climate Variation
ggplot(sister_species_coords, aes(x = temperature, y = precipitation)) +
  geom_point(aes(color = species_name), alpha = 0.6, size = 3) +
  labs(
    title = "Climate Variability Across Sister Species Locations",
    x = "Annual Mean Temperature (°C)",
    y = "Annual Precipitation (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_viridis_d(option = "turbo") # Adds color differentiation for species

# Interpretation: This plot helps investigate whether climate differences correlate with geographic speciation patterns.
