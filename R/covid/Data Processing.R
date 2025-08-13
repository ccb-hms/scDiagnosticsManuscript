# -----------------------------------------------------
# COVID-19 - Data Processing
# -----------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scran)
library(scater)

# Source files
source("R/auxiliary/addReferencePCA.R")

# -----------------------------------------------------

# Load the formatted SCE objects
normal_data <- readRDS("data/covid/normal_data.rds")
covid_data <- readRDS("data/covid/covid_data.rds")

# Setting the seed
set.seed(0)

# Remove disease column from single-disease datasets
colData(normal_data)$disease <- NULL
colData(covid_data)$disease <- NULL

# Apply PCA to reference SCE object using union HVGs (reference & query)
normal_data <- addReferencePCA(normal_data, covid_data, "normal")

# For COVID data, use its own HVGs (since it won't be used as reference)
covid_dec <- suppressWarnings(modelGeneVar(covid_data))
covid_hvgs <- getTopHVGs(covid_dec, n = 2000)
covid_data <- runPCA(covid_data, subset_row = covid_hvgs, ncomponents = 50)

# Set the annotation data as character data
normal_data$cell_type <- as.character(normal_data$cell_type)
covid_data$cell_type <- as.character(covid_data$cell_type)

# Save all processed datasets
saveRDS(normal_data, "data/covid/normal_data.rds")
saveRDS(covid_data, "data/covid/covid_data.rds") 

cat("Processing complete! Saved:\n")
cat("- normal_data.rds:", ncol(normal_data), "cells\n")
cat("- covid_data.rds:", ncol(covid_data), "cells\n") 

# Clean up
rm(normal_data, covid_data)
