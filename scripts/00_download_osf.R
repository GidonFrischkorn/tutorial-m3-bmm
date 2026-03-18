# ============================================================================
# Download fitted model objects, figures, and data from OSF
# ============================================================================

pacman::p_load(osfr, here)
library(osfr)
library(here)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
osf_project_id <- "Yb7wm"

# ---------------------------------------------------------------------------
# Retrieve public project
# ---------------------------------------------------------------------------
message("Connecting to OSF project: ", osf_project_id)
project <- osf_retrieve_node(osf_project_id)
message("Project: ", project$name, "\n")

# ---------------------------------------------------------------------------
# Helper: list all files in an OSF folder
# ---------------------------------------------------------------------------
get_osf_folder <- function(project, folder_name) {
  folders <- osf_ls_files(project, type = "folder")
  match <- folders[folders$name == folder_name, ]
  if (nrow(match) == 0) {
    warning("Folder '", folder_name, "' not found on OSF. Skipping.")
    return(NULL)
  }
  match
}

# ---------------------------------------------------------------------------
# Download output/ files
# ---------------------------------------------------------------------------
osf_output <- get_osf_folder(project, "output")

if (!is.null(osf_output)) {
  output_dir <- here("output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_files <- osf_ls_files(osf_output)
  n_files <- nrow(output_files)

  message("--- Downloading fitted model objects (output/) ---")
  message("  Found ", n_files, " files on OSF\n")

  for (i in seq_len(n_files)) {
    fname <- output_files$name[i]
    local_path <- file.path(output_dir, fname)

    if (file.exists(local_path)) {
      message("  SKIP (already exists): ", fname)
      next
    }

    message("  Downloading: ", fname, " (", i, "/", n_files, ")")
    osf_download(output_files[i, ], path = output_dir, conflicts = "skip")
  }
}

# ---------------------------------------------------------------------------
# Download data/ files
# ---------------------------------------------------------------------------
osf_data <- get_osf_folder(project, "data")

if (!is.null(osf_data)) {
  data_dir <- here("data")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  data_files <- osf_ls_files(osf_data)
  n_files <- nrow(data_files)

  message("\n--- Downloading datasets (data/) ---")
  message("  Found ", n_files, " files on OSF\n")

  for (i in seq_len(n_files)) {
    fname <- data_files$name[i]
    local_path <- file.path(data_dir, fname)

    if (file.exists(local_path)) {
      message("  SKIP (already exists): ", fname)
      next
    }

    message("  Downloading: ", fname, " (", i, "/", n_files, ")")
    osf_download(data_files[i, ], path = data_dir, conflicts = "skip")
  }
}

# ---------------------------------------------------------------------------
# Download figures/ files
# ---------------------------------------------------------------------------
osf_figures <- get_osf_folder(project, "figures")

if (!is.null(osf_figures)) {
  figures_dir <- here("figures")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  figures_files <- osf_ls_files(osf_figures)
  n_files <- nrow(figures_files)

  message("\n--- Downloading figures (figures/) ---")
  message("  Found ", n_files, " files on OSF\n")

  for (i in seq_len(n_files)) {
    fname <- figures_files$name[i]
    local_path <- file.path(figures_dir, fname)

    if (file.exists(local_path)) {
      message("  SKIP (already exists): ", fname)
      next
    }

    message("  Downloading: ", fname, " (", i, "/", n_files, ")")
    osf_download(figures_files[i, ], path = figures_dir, conflicts = "skip")
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
output_files_local <- list.files(here("output"), pattern = "\\.rds$",
                                 full.names = TRUE)
figures_local <- list.files(here("figures"), pattern = "\\.pdf$",
                            full.names = TRUE)
total_mb <- round(sum(file.size(c(output_files_local, figures_local))) / 1e6, 1)

message("\nDone. Local output/ contains ", length(output_files_local),
        " .rds files; figures/ contains ", length(figures_local), " .pdf files",
        " (", total_mb, " MB total).")
message("You can now render the manuscript without re-fitting models.")
