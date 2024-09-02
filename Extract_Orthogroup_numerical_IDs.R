# Load necessary libraries
library(dplyr)
library(readr)

# List of GFF files (one for each species)
gff_files <- list(
  "Lmajor_Standardised.gff",
  "Angomonas_Standardised.gff",
  "Crithidia_Standardised.gff",
  "Lot6_Standardised.gff",
  "Porcisia_Standardised.gff",
  "Tbrucei4_Standardised.gff",
  "Tbrucei8_Standardised.gff",
  "Tcruzi_Standardised.gff",
  "Tryp_Standardised.gff"
)

# Step 1: Create a global mapping from orthogroup IDs to numerical values
# Combine orthogroup IDs from all GFF files to ensure a consistent global mapping
all_orthogroups <- data.frame()

for (gff_file in gff_files) {
  gff_data <- read_tsv(gff_file, col_names = FALSE)
  colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description", "gene_name", "Orthogroup")
  all_orthogroups <- bind_rows(all_orthogroups, gff_data %>% select(Orthogroup))
}

# Step 2: Extract and sort the orthogroup IDs
unique_orthogroups <- unique(all_orthogroups$Orthogroup)
sorted_orthogroups <- sort(unique_orthogroups, na.last = TRUE)  # Sort orthogroups and keep NA at the end

# Step 3: Create a global mapping from sorted orthogroup IDs to numerical values
orthogroup_map <- setNames(seq_along(sorted_orthogroups), sorted_orthogroups)

# Save the mapping to a file for reference (optional)
write_csv(data.frame(Orthogroup = sorted_orthogroups, NumericalID = seq_along(sorted_orthogroups)), "orthogroup_mapping.csv")

# Step 4: Process each updated GFF file to extract numerical orthogroup IDs
process_gff_file <- function(gff_file, orthogroup_map) {
  gff_data <- read_tsv(gff_file, col_names = FALSE)
  colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description", "gene_name", "Orthogroup")
  
  # Convert orthogroup IDs to numerical values using the global mapping
  gff_data <- gff_data %>%
    mutate(num_orthogroup = orthogroup_map[Orthogroup])
  
  # Extract the numerical orthogroup IDs in the exact order and write to a text file
  numerical_orthogroup_ids <- gff_data %>%
    select(num_orthogroup) %>%
    filter(!is.na(num_orthogroup)) %>%
    pull()
  
  # Define the output file name based on the input GFF file name
  output_file <- sub("\\.gff$", "_numerical_orthogroups.txt", gff_file)
  write(paste(numerical_orthogroup_ids, collapse = " "), file = output_file)
}

# Apply the function to each updated GFF file
for (gff_file in gff_files) {
  process_gff_file(gff_file, orthogroup_map)
}



