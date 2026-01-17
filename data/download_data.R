# data/download_data.R
# Usage: source("data/download_data.R"); download_all_data()

download_all_data <- function() {
  
  # Zenodo dataset: https://doi.org/10.5281/zenodo.18274942
  zenodo_base <- "https://zenodo.org/records/18274942/files/"
  
  zenodo_urls <- list(
    covid_sce = paste0(zenodo_base, "covid_data_sce.rds"),
    covid_ref = paste0(zenodo_base, "normal_data_sce.rds"),
    merfish_sce = paste0(zenodo_base, "dss9_data.rds"),
    merfish_ref = paste0(zenodo_base, "healthy_data.rds")
  )
  
  # Create subdirectories if needed
  dir.create("data/covid", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/merfish", showWarnings = FALSE, recursive = TRUE)
  
  # Download COVID data
  message("Downloading COVID-19 data...")
  download.file(zenodo_urls$covid_sce, "data/covid/covid_data_sce.rds", mode = "wb")
  download.file(zenodo_urls$covid_ref, "data/covid/normal_data_sce.rds", mode = "wb")
  
  # Download MERFISH data
  message("Downloading MERFISH data...")
  download.file(zenodo_urls$merfish_sce, "data/merfish/dss9_data.rds", mode = "wb")
  download.file(zenodo_urls$merfish_ref, "data/merfish/healthy_data.rds", mode = "wb")
  
  message("✓ All data downloaded successfully!")
  message("Data saved to data/covid/ and data/merfish/")
}

# Run this to download
download_all_data()