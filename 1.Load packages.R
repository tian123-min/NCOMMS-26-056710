# Load and install required packages.

candidate_libraries <- c(
  Sys.getenv("NCOMMS_R_LIB", unset = NA_character_),
  file.path(Sys.getenv("OUTPUT_DIR", "outputs"), "package_library"),
  file.path(tempdir(), "ncomms_r_library")
)
candidate_libraries <- candidate_libraries[!is.na(candidate_libraries) & nzchar(candidate_libraries)]

local_library <- NULL
for (candidate in candidate_libraries) {
  dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(candidate) && file.access(candidate, 2) == 0) {
    local_library <- candidate
    break
  }
}

if (is.null(local_library)) {
  stop("No writable R package library is available. Set NCOMMS_R_LIB to a writable directory.")
}

.libPaths(c(local_library, .libPaths()))

required_packages <- c(
  "ranger",
  "caret",
  "tidyverse",
  "readxl",
  "writexl",
  "openxlsx",
  "matrixStats",
  "h2o",
  "e1071",
  "RSNNS",
  "doParallel"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, lib = local_library, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

theme_set(theme_bw())
