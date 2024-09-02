# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

## This code intends to reverse-map the inferred ancestor output from MLGO from numerical ID to OG IDs.
## This code requires the 'orthogroup_mapping.csv' created in Code3, which was used for inference for accurate reverse-mapping.

# Load the orthogroup mapping CSV
orthogroup_mapping <- read.csv("orthogroup_mapping.csv")

# Create a lookup table for mapping NumericalID to Orthogroup
id_to_orthogroup <- setNames(orthogroup_mapping$Orthogroup, orthogroup_mapping$NumericalID)

# Read the ancestral gene order.out file
geneorder_file <- "geneorder.out"
geneorder_lines <- readLines(geneorder_file)

# Initialize an empty list to store gene orders
ancestral_gene_orders <- list()

# Parse the file and populate the list
current_key <- NULL
for (line in geneorder_lines) {
  if (startsWith(line, ">")) {
    # Extract the key (name) for the gene order
    current_key <- sub(">", "", line)
    ancestral_gene_orders[[current_key]] <- c()
  } else if (!is.null(current_key)) {
    # Append the numerical gene data to the current list entry
    gene_numbers <- as.numeric(strsplit(line, " ")[[1]])
    ancestral_gene_orders[[current_key]] <- c(ancestral_gene_orders[[current_key]], gene_numbers)
  }
}

# Function to map numerical IDs to orthogroups and orientations
map_genes <- function(ancestral_genes, id_to_orthogroup) {
  sapply(abs(ancestral_genes), function(gene_id) {
    orthogroup <- id_to_orthogroup[as.character(gene_id)]
    orthogroup
  })
}

# Process each ancestral gene order and save to CSV
for (name in names(ancestral_gene_orders)) {
  ancestral_genes <- ancestral_gene_orders[[name]]
  mapped_genes <- map_genes(ancestral_genes, id_to_orthogroup)
  gene_order_df <- data.frame(order = seq_along(mapped_genes), gene = mapped_genes)
  output_file <- paste0("ancestral_gene_order_", name, ".csv")
  write.csv(gene_order_df, output_file, row.names = FALSE)
  message("Saved: ", output_file)
}
