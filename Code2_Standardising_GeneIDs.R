# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)

## Standardising Angomonas specific gene ID to respective OG ID

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr2_Angomonas.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*Name=([^;]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Angomonas_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Angomonas_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising Tbrucei_Chr4 specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr4_Tbrucei_927_editted.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^:]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Tbrucei4_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Tbrucei4_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising Tbrucei_Chr8 specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr8_Tbrucei_927_editted.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^:]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Tbrucei8_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Tbrucei8_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)



## Standardising Lotmaria_Chr5 specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr5_Lotmaria.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^;]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Lotmaria5_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Lot5_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)



## Standardising Lotmaria_Chr6 specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr6_Lotmaria.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^;]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Lotmaria6_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Lot6_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising Porcisia specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr31_Porcisia.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^-CDS]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Porcisia_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Porcisia_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising T.cruzi specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr31_Tcruz_NONEsmeraldo.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^:]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Tcruzi_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Tcruzi_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising Crithidia specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr35_Crithidia.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^;]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Crithidia_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Crithidia_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising Tryp_sp specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Contig15_TrypSp.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*ID=([^;]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Tryp_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Tryp_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)


## Standardising L.major specific gene ID to respective OG ID 

# Step 1: Read the GFF file and extract the gene name
gff_data <- read_tsv("Chr31_LMajor.gff", comment = "#", col_names = FALSE)
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "description")

# Extract gene name from the attributes column
gff_data <- gff_data %>%
  mutate(gene_name = sub(".*gene-([^:]+).*", "\\1", description))

# Step 2: Read the orthogroup.tsv file
orthogroup_file <- "Lmajor_Orthogroup.tsv"
orthogroup_data <- read_tsv(orthogroup_file, col_names = FALSE)

# Convert orthogroup data to long format
orthogroup_long <- orthogroup_data %>%
  pivot_longer(cols = -X1, values_to = "gene_name") %>%
  filter(!is.na(gene_name)) %>%
  rename(Orthogroup = X1)

# Step 3: Match gene names to orthogroup ID and add a new column in the GFF data
gff_data <- gff_data %>%
  left_join(orthogroup_long, by = "gene_name", keep = TRUE)

# Ensure original order is maintained
gff_data <- gff_data %>%
  arrange(row_number())

# Step 4: Write the updated GFF data to a new file
output_file <- "Lmajor_Standardised.gff"
write_tsv(gff_data, output_file, col_names = FALSE)
