# data/download_data.R
# Usage: source("data/download_data.R")

download_all_data <- function() {
  
  zenodo_urls <- list(
    covid_sce = "https://zenodo.org/.../covid_data_sce.rds",
    covid_ref = "https://zenodo.org/.../normal_data_sce.rds",
    merfish_sce = "https://zenodo.org/.../dss9_data.rds",
    merfish_ref = "https://zenodo.org/.../healthy_data.rds"
  )
  
  # Create subdirectories if needed
  dir.create("data/covid", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/merfish", showWarnings = FALSE, recursive = TRUE)
  
  # Download COVID data
  message("Downloading COVID-19 data...")
  download.file(zenodo_urls$covid_sce, "data/covid/covid_data_sce.rds")
  download.file(zenodo_urls$covid_ref, "data/covid/normal_data_sce.rds")
  
  # Download MERFISH data
  message("Downloading MERFISH data...")
  download.file(zenodo_urls$merfish_sce, "data/merfish/dss9_data.rds")
  download.file(zenodo_urls$merfish_ref, "data/merfish/healthy_data.rds")
  
  message("✓ All data downloaded successfully!")
}