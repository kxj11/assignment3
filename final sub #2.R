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
library(dplyr)
library(readr) # For read_tsv()

# Data Loading and Initial Filtering -------------------------------------------
# Load and filter data (removing rows with missing values) from BOLD dataset
# Load data and convert to data frame 
#dflizard <- read_tsv("C:/Users/Dell User/Documents/BINF/bold_data-2.tsv", show_col_types = FALSE) %>%
#  as.data.frame()
#experimen:
setwd("../6210 assignment 2")

dflizard <- read_tsv("../6210 assignment 2/bold_data-2.tsv", show_col_types = FALSE) %>%
  as.data.frame()


# This dataset includes species names, countries, latitude, longitude, and COI sequences.
# Now apply the select and filter functions on the data frame
dflizardfilter <- dflizard %>%
  dplyr::select(processid, species_name, country, lat, lon, nucleotides) %>%
  dplyr::filter(!is.na(species_name), !is.na(lat), !is.na(lon), !is.na(nucleotides), !is.na(country))

# View the filtered data to confirm
print(dflizardfilter)

# Fetch COI Gene Sequences from NCBI -------------------------------------------
# Fetch sequence data from NCBI for Phelsuma species to build a phylogenetic tree.
# Saves sequences in FASTA format for reproducibility and future use.
# Search for Phelsuma COI gene sequences
search <- entrez_search(db = "nucleotide", term = "Phelsuma COI", retmax = 100)

# Proceed if there are results
if (search$count > 0) {
  
  # Retrieve metadata for the fetched sequence IDs
  metadata <- entrez_summary(db = "nucleotide", id = search$ids)
  
  # Extract organism names from the metadata
  species_info <- sapply(metadata, function(x) x$organism)  # Use the 'organism' field for species names
  
  # Create a data frame with species names
  sequence_data <- data.frame(Species = species_info, stringsAsFactors = FALSE)
  
  # View the combined data
  print(sequence_data)
  
  # saving them as Fasta files 
  
  # Retrieve the COI sequences in FASTA format using sequence IDs
  fasta_sequences <- entrez_fetch(db = "nucleotide", id = search$ids, rettype = "fasta", retmode = "text")
  
  # Save the sequences to a FASTA file on your computer
  write(fasta_sequences, file = "../6210 assignment 2/phelsuma_sequences.fasta")
  
  # View the first few lines of the saved FASTA file (optional)
  cat(readLines("../6210 assignment 2/phelsuma_sequences.fasta", n = 10), sep = "\n")
  
  cat("FASTA sequences saved successfully!")
}

# Phylogenetic Analysis --------------------------------------------------------
# Load sequence data, align, and construct phylogenetic tree for species comparison.
# This helps visualize evolutionary relationships among species and test speciation patterns.
#need to be aligned 
install.packages("BiocManager")
BiocManager::install("msa")
library(msa)
alignment <- msa::msa("../6210 assignment 2/phelsuma_sequences.fasta", type = "dna")
print(alignment)

# Read the alignment directly
alignment <- msa::msa("../6210 assignment 2/phelsuma_sequences.fasta", type = "dna")

# If needed, convert the alignment to DNAbin format (check if it's already DNAbin first)
if (!inherits(alignment, "DNAbin")) {
  alignment_dnabin <- as.DNAbin(alignment)
} else {
  alignment_dnabin <- alignment
}

# Compute the distance matrix and construct the phylogenetic tree using Neighbor-Joining
dist_matrix <- dist.dna(alignment_dnabin)
tree <- nj(dist_matrix)

# Assuming `sequence_data` has the species names in a column named 'Species'
species_names <- sequence_data$Species
tree$tip.label <- species_names

# Set plot margins to ensure the title is fully visible
par(mar = c(5, 1, 9, 1))  # Further increased top margin (bottom, left, top, right)

# Plot the tree with improved readability and full title display
plot(
  tree,
  main = "Phylogenetic Tree of Phelsuma Species",
  cex = 0.5,             # Smaller label size
  label.offset = 0.02,   # Space between tips and labels
  no.margin = TRUE       # Removes outer margins to use full plot area
)

# Interpretation: This tree provides insight into which species are closely related (sister species),
# which addresses the main question of speciation patterns in *Phelsuma*

# Combine BOLD and NCBI Data for Geographic Analysis ---------------------------
# Check structure of both datasets
str(dflizardfilter)
str(sequence_data)

# Perform the left join with specified column names
combined_data <- dflizardfilter %>%
  left_join(sequence_data, by = c("species_name" = "Species"))

#it is found to have duplicates, So we will aggregate the dflizardfilter data

# dflizardfilter is already defined
# Aggregate data by species_name and country
aggregated_data <- dflizardfilter %>%
  group_by(species_name, country) %>%
  summarise(
    mean_lat = mean(lat, na.rm = TRUE),     # Calculate mean latitude
    mean_lon = mean(lon, na.rm = TRUE),     # Calculate mean longitude
    sequence_count = n()                     # Count number of sequences
  )

# Now perform the left join with sequence_data
combined_data <- aggregated_data %>%
  left_join(sequence_data, by = c("species_name" = "Species"))

# View the combined data
print(combined_data)

# still finding duplicates in aggregate data and sequence data so we are checking for duplicates again

# Check for duplicates in aggregated_data
duplicates_aggregated <- aggregated_data %>%
  group_by(species_name) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# Check for duplicates in sequence_data
duplicates_sequence <- sequence_data %>%
  group_by(Species) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# Print results
print("Duplicates in Aggregated Data:")
print(duplicates_aggregated)

print("Duplicates in Sequence Data:")
print(duplicates_sequence)

#checked for duplicates above and now rectifying that below:

# Aggregate dflizardfilter to remove duplicates by calculating means
aggregated_data <- dflizardfilter %>%
  group_by(species_name, country) %>%
  summarise(
    mean_lat = mean(lat, na.rm = TRUE),
    mean_lon = mean(lon, na.rm = TRUE),
    sequence_count = n(),
    .groups = "drop"  # This removes the grouping after summarizing
  )

# Aggregate sequence_data to remove duplicates by keeping unique species
unique_sequence_data <- sequence_data %>%
  group_by(Species) %>%
  summarise(count = n(), .groups = "drop")  # Count occurrences

# Join the aggregated datasets
combined_data <- aggregated_data %>%
  left_join(unique_sequence_data, by = c("species_name" = "Species"))

# View the combined data
print(combined_data)

# Phylogenetic Tree Cleanup ----------------------------------------------------
# Remove duplicate tips and clean up tree for better visualization.
#now remove duplicates further

# Get unique tip labels
unique_labels <- unique(tree$tip.label)

# Create a new tree with unique labels
# Assuming you want to keep the first occurrence
tree_unique <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique_labels)])

# Ensure that the underlying data used to create the tree does not have duplicates.
# Check the structure of the new tree
print(tree_unique)

# Plot the tree to ensure it looks correct
plot(tree_unique, main = "Phylogenetic Tree of Phelsuma Species (Unique Labels)")
#making it look better
plot(tree_unique, cex = 0.5, main = "Phylogenetic Tree of Phelsuma Species (Unique Labels)")

# Sister Species Identification ------------------------------------------------
# Identify sister species pairs in the phylogenetic tree.
# Step 1: Create a function to identify sister species
find_sister_species <- function(tree) {
  sister_pairs <- list()
  
  # Loop over each internal node in the tree
  for (i in (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))) {
    # Get the descendants of the current node
    descendants <- tree$tip.label[Descendants(tree, i, type = "tips")[[1]]]
    
    # If there are exactly two descendants, we have a sister species pair
    if (length(descendants) == 2) {
      sister_pairs[[length(sister_pairs) + 1]] <- descendants
    }
  }
  
  # Return the sister pairs as a data frame
  sister_pairs_df <- do.call(rbind, sister_pairs) %>% as.data.frame()
  colnames(sister_pairs_df) <- c("species_1", "species_2")  # Name the columns
  
  return(sister_pairs_df)
}

# Step 2: Identify sister species
sister_pairs_df <- find_sister_species(tree_unique)

print(sister_pairs_df)
class(combined_data)
colnames(combined_data)

#filtering for self-pairs and making it clean
sister_pairs_df <- sister_pairs_df %>%
  filter(species_1 != species_2)

# Check for self-pairs directly
self_pairs <- sister_pairs_df %>% filter(species_1 == species_2)
self_pairs <- sister_pairs_df %>% filter(species_1 == species_2)
if (nrow(self_pairs) == 0) {
  cat("No self-pairs found.\n")
} else {
  print(self_pairs)
}

# Step 3: Extract geographic data for sister species
geographic_sisters <- combined_data %>%
  filter(species_name %in% unlist(c(sister_pairs_df$species_1, sister_pairs_df$species_2))) %>%
  dplyr::select(species_name, country, mean_lat, mean_lon)

# View geographic data for sister species
print(geographic_sisters)

# Plot geographic distribution of sister species to visualize spatial patterns in speciation.
# Step 4: Visualization
ggplot(data = geographic_sisters, aes(x = mean_lon, y = mean_lat, label = species_name)) +
  geom_point(aes(color = country, shape = country), size = 3) +  # Different shapes and colors for each country
  geom_text_repel(
    size = 2,                         # Smaller font size for labels
    fontface = "bold",                # Bold font for readability
    box.padding = 0.5,                # Padding around the labels
    point.padding = 0.3,              # Space between points and labels
    segment.color = "grey",           # Connecting lines to points
    segment.size = 0.5,               # Thickness of connecting lines
    max.overlaps = 50                 # Allow more overlapping labels to be shown
  ) +
  borders("world") +
  labs(
    title = "Geographic Distribution of Sister Species in Phelsuma",
    x = "Longitude",
    y = "Latitude",
    color = "Country",
    shape = "Country"
  ) +
  theme_minimal()

# Interpretation: This plot helps reveal whether sister species occupy the same or different regions,
# providing insights into potential allopatric or sympatric speciation.

# answering the research question further.
# Categorize Sister Species by Region Overlap--------------------------------------
# Aggregate `geographic_sisters` to remove duplicate species entries by calculating average coordinates
geographic_sisters_unique <- geographic_sisters %>%
  group_by(species_name) %>%
  summarise(
    mean_lat = mean(mean_lat, na.rm = TRUE),
    mean_lon = mean(mean_lon, na.rm = TRUE),
    country = first(country)  # Assumes a single country per species
  ) %>%
  ungroup()

# Prepare `sister_lines_data` by merging unique sister species with their geographic coordinates
sister_lines_data <- sister_pairs_df %>%
  left_join(geographic_sisters_unique, by = c("species_1" = "species_name")) %>%
  rename(lat_1 = mean_lat, lon_1 = mean_lon, country_1 = country) %>%
  left_join(geographic_sisters_unique, by = c("species_2" = "species_name")) %>%
  rename(lat_2 = mean_lat, lon_2 = mean_lon, country_2 = country)

# Filter out rows with missing geographic information in `sister_lines_data`
sister_lines_data_clean <- sister_lines_data %>%
  filter(!is.na(lat_1) & !is.na(lon_1) & !is.na(lat_2) & !is.na(lon_2))

# View the cleaned data to confirm missing entries are removed
head(sister_lines_data_clean)


# Create a summary table for sister species pairs with region labels
sister_summary_table <- sister_lines_data_clean %>%
  mutate(region_overlap = ifelse(country_1 == country_2, "Same Region", "Different Region")) %>%
  dplyr::select(species_1, species_2, country_1, country_2, region_overlap)

# View the summary table
print(sister_summary_table)

#cleaning the above table
# Filter to show only pairs with complete region data
sister_summary_table_complete <- sister_summary_table %>%
  filter(!is.na(region_overlap))

# Print the filtered table
print(sister_summary_table_complete)

#making the table visually better
# Install knitr
# Display the table in a clean format
kable(sister_summary_table_complete, 
      col.names = c("Species 1", "Species 2", "Country 1", "Country 2", "Region Overlap"),
      caption = "Geographic Overlap of Sister Species in Phelsuma",
      align = "l")  # left-align columns

#Quantify Region Overlap--------------------------------------------------------------
# Summarize and categorize sister species based on geographic overlap to quantify speciation patterns.
# Summary of region overlap counts
overlap_summary <- sister_summary_table_complete %>%
  group_by(region_overlap) %>%
  summarize(count = n())

# Print the summary
print(overlap_summary)

# make it better looking and adding a visualization
overlap_summary <- as.data.frame(overlap_summary)
colnames(overlap_summary) <- c("Geographic Overlap", "Number of Pairs")

print(overlap_summary)

#check structure
str(overlap_summary)

# Interpretation: The table quantifies the occurrence of speciation in shared vs. different regions,
# addressing the main question of geographic influence on speciation.

##calculate geographical distances-------------------------------------------------------
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
options(repr.plot.width=14, repr.plot.height=9)

# Generate the improved plot
ggplot(data = sister_lines_data_clean, aes(x = lon_1, y = lat_1)) +
  borders("world", colour = "grey85", fill = "grey95") +
  geom_segment(aes(xend = lon_2, yend = lat_2, color = distance_category, 
                   linetype = distance_category, size = distance_category), alpha = 0.5) +  # Reduced alpha for lines
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
    max.overlaps = Inf, segment.size = 0.3, box.padding = 0.4, force = 1  # Adjusted for readability
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

#check patterns in my data----------------------------------------------------------
#1. Test if Sister Species are More Likely in the Same or Different Regions
# Count table for "Same Region" and "Different Region"
region_table <- table(sister_lines_data_clean$distance_category)

# Run the Chi-square test
chi_test <- chisq.test(region_table)

# View the test results
print(chi_test)

# Interpretation: The Chi-square test evaluates whether sister species pairs are more likely
# to exist in the same or different regions, addressing the main hypothesis on speciation.The p-value (0.2) was insignificant suggesting the observed distribution of sister species pairs across regions could be due to chance rather than a clear preference for same or different regions

#adding and exploring using the climate data-----------------------------------------------------------------

# Define the path where the .tif file is located after extraction
bio1_path <- "../6210 assignment 2/wc2.1_2.5m_bio_1/wc2.1_2.5m_bio_1.tif"

# Load the raster
bio1 <- raster(bio1_path)

#precipitation info
# Define the path where the .tif file is located after extraction
bio12_path <- "../6210 assignment 2/wc2.1_2.5m_bio_12/wc2.1_2.5m_bio_12.tif"

# Load the raster
bio12 <- raster(bio12_path)

#climate analysis-------------------------------------------------------------------------
#Code to Create sister_species_coords
# Get a list of unique sister species names from `sister_pairs_df`
sister_species <- unique(c(sister_pairs_df$species_1, sister_pairs_df$species_2))

#  Filter `dflizardfilter` to include only sister species
sister_species_coords <- dflizardfilter %>%
  filter(species_name %in% sister_species) %>%
  dplyr::select(species_name, lat, lon) %>%
  dplyr::rename(latitude = lat, longitude = lon)

# View the `sister_species_coords` data frame to confirm
print(sister_species_coords)

#visualization
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
