# ---------------------------------
# Environment Setup for CellTypist
# ---------------------------------

# Load libraries
library(reticulate)

# ---------------------------------

# Set up Python environment for CellTypist
environmentSetupCellTypist <- function(){

  cat("=== Python Environment Setup ===\n")

  # Function to safely print Python info
  print_python_info <- function() {
    tryCatch({
      # Method 1: Try py_config
      config <- py_config()
      if(!is.null(config)) {
        cat("Python executable:", 
            if(is.character(config$python)) config$python else "Found", "\n")
      }
    }, error = function(e) {
      # Method 2: Try direct Python query
      tryCatch({
        py_run_string("
  import sys
  version_info = '.'.join(map(str, sys.version_info[:3]))
  print('Python version:', version_info)
  print('Python path:', sys.executable)
  ")
      }, error = function(e2) {
        cat("Python information not available\n")
      })
    })
  }

  # Try to set up conda environment
  setup_conda_env <- function() {
    cat("Setting up conda environment...\n")
    
    # Check for conda
    has_conda <- FALSE
    tryCatch({
      conda_path <- conda_binary()
      if(!is.null(conda_path) && conda_path != "" && file.exists(conda_path)) {
        has_conda <- TRUE
        cat("✓ Found conda\n")
      }
    }, error = function(e) {
      has_conda <- FALSE
    })
    
    # Install miniconda if needed
    if(!has_conda) {
      cat("Installing miniconda...\n")
      tryCatch({
        install_miniconda()
        has_conda <- TRUE
        cat("✓ Miniconda installed\n")
        Sys.sleep(2)  # Give it time to set up
      }, error = function(e) {
        cat("✗ Miniconda install failed\n")
        return(FALSE)
      })
    }
    
    if(has_conda) {
      env_name <- "celltypist_env"
      
      # Create environment if needed
      tryCatch({
        existing_envs <- conda_list()$name
        if(!env_name %in% existing_envs) {
          cat("Creating environment:", env_name, "\n")
          conda_create(env_name, python_version = "3.9")
          
          # Install packages
          cat("Installing packages...\n")
          conda_install(env_name, c("pandas", "numpy", "scipy", "h5py"))
          conda_install(env_name, "scanpy", channel = "conda-forge")
          conda_install(env_name, "celltypist", pip = TRUE)
          cat("✓ Packages installed\n")
        }
        
        # Use environment
        use_condaenv(env_name, required = TRUE)
        return(TRUE)
        
      }, error = function(e) {
        cat("✗ Environment setup failed:", e$message, "\n")
        return(FALSE)
      })
    }
    
    return(FALSE)
  }

  # Main setup
  success <- setup_conda_env()

  if(!success) {
    cat("Trying system Python...\n")
    tryCatch({
      py_install(c("pandas", "numpy", "scanpy", "celltypist"))
      success <- TRUE
      cat("✓ System Python setup complete\n")
    }, error = function(e) {
      cat("✗ System Python setup failed\n")
    })
  }

  # Print Python info
  print_python_info()

  # Test imports
  cat("\nTesting imports...\n")
  test_packages <- c("pandas", "numpy", "scanpy", "celltypist")
  imported_packages <- character(0)

  for(pkg in test_packages) {
    tryCatch({
      import(pkg)
      cat("✓", pkg, "\n")
      imported_packages <- c(imported_packages, pkg)
    }, error = function(e) {
      cat("✗", pkg, ":", e$message, "\n")
    })
  }

  # Summary
  cat("\n=== Setup Summary ===\n")
  cat("Successfully imported:", length(imported_packages), "/", length(test_packages), "packages\n")

  if("celltypist" %in% imported_packages && "scanpy" %in% imported_packages) {
    cat("🎉 Ready for CellTypist pipeline!\n")
  } else {
    cat("⚠️  Setup incomplete. You may need to install packages manually.\n")
    cat("\nManual installation commands:\n")
    cat("conda create -n celltypist_env python=3.9\n")
    cat("conda activate celltypist_env\n") 
    cat("conda install scanpy pandas numpy -c conda-forge\n")
    cat("pip install celltypist\n")
  }

}

