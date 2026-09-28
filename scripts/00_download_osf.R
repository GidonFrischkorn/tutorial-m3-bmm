# Downloads the fitted model objects (output/), data (data/), and figures
# (figures/) from the OSF project of this tutorial, so the tutorial scripts and
# the manuscript run without refitting. Existing local files are skipped.
# Sourced optionally at the top of each tutorial script.

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(osfr, here)

osf_project_id <- "Yb7wm"

message("Connecting to OSF project: ", osf_project_id)
project <- osf_retrieve_node(osf_project_id)
message("Project: ", project$name, "\n")

# Returns the OSF folder with the given name, or NULL if it does not exist.
# n_max = Inf: osf_ls_files() lists only 10 entries by default.
get_osf_folder <- function(project, folder_name) {
  folders <- osf_ls_files(project, type = "folder", n_max = Inf)
  match <- folders[folders$name == folder_name, ]
  if (nrow(match) == 0) {
    warning("Folder '", folder_name, "' not found on OSF. Skipping.")
    return(NULL)
  }
  match
}

###############################################################################!
# 1) Download Files ------------------------------------------------------------
###############################################################################!

## 1.1) Fitted model objects (output/) -----------------------------------------
osf_output <- get_osf_folder(project, "output")

if (!is.null(osf_output)) {
  output_dir <- here("output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_files <- osf_ls_files(osf_output, n_max = Inf)
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

## 1.2) Data (data/) -----------------------------------------------------------
osf_data <- get_osf_folder(project, "data")

if (!is.null(osf_data)) {
  data_dir <- here("data")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  data_files <- osf_ls_files(osf_data, n_max = Inf)
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

## 1.3) Figures (figures/) -----------------------------------------------------
osf_figures <- get_osf_folder(project, "figures")

if (!is.null(osf_figures)) {
  figures_dir <- here("figures")
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  figures_files <- osf_ls_files(osf_figures, n_max = Inf)
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

###############################################################################!
# 2) Summary -------------------------------------------------------------------
###############################################################################!

local_dirs  <- c("output", "data", "figures")
local_files <- lapply(local_dirs, function(d) {
  list.files(here(d), full.names = TRUE)
})
total_mb <- round(sum(file.size(unlist(local_files))) / 1e6, 1)

message("\nDone. Local files: ",
        paste0(lengths(local_files), " in ", local_dirs, "/", collapse = ", "),
        " (", total_mb, " MB total).")
