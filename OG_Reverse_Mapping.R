# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the orthogroup mapping CSV
orthogroup_mapping <- read.csv("orthogroup_mapping.csv")

# Create a lookup table for mapping NumericalID to Orthogroup
id_to_orthogroup <- setNames(orthogroup_mapping$Orthogroup, orthogroup_mapping$NumericalID)

# Read the ancestral gene order.out file
geneorder_file <- "run3_geneorder.out"
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
  output_file <- paste0("Run3_ancestral_gene_order_", name, ".csv")
  write.csv(gene_order_df, output_file, row.names = FALSE)
  message("Saved: ", output_file)
}

# Print the first few rows of one of the output files to verify
head(gene_order_df)

____

# Assign colours to each orthogroup
orthogroup_colours <- setNames(rainbow(length(unique(orthogroup_mapping$Orthogroup))), unique(orthogroup_mapping$Orthogroup))

# Create a function to plot gene orders
plot_gene_order <- function(ancestral_genes, name, id_to_orthogroup, orthogroup_colours) {
  mapped_genes <- map_genes(ancestral_genes, id_to_orthogroup)
  gene_order_df <- data.frame(order = seq_along(mapped_genes), gene = mapped_genes)
  
  gene_order_df$Orthogroup <- gsub("^[+-]", "", gene_order_df$gene)
  gene_order_df$Colour <- orthogroup_colours[gene_order_df$Orthogroup]
  
  ggplot(gene_order_df, aes(x = order, y = 1, fill = Colour)) +
    geom_tile() +
    scale_fill_identity() +
    labs(title = paste("Ancestral Gene Order:", name), x = "Gene Order", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

# Plot and save each ancestral gene order
for (name in names(ancestral_gene_orders)) {
  ancestral_genes <- ancestral_gene_orders[[name]]
  plot <- plot_gene_order(ancestral_genes, name, id_to_orthogroup, orthogroup_colours)
  ggsave(paste0("ancestral_gene_order_", name, ".png"), plot)
  print(plot)
}
