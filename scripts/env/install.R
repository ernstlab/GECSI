#!/usr/bin/env Rscript
set -e

rm(list = ls(all = TRUE))
options(warn = -1)
options(error = function() {
  traceback(2)
  if (!interactive()) quit("no", status = 1, runLast = FALSE)
})
gc()
if (.Platform$OS.type == "windows") memory.limit(size = Inf)
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")


library(readr)
library(pacman)

env_file <- "./scripts/env/R_environment.txt"
if (!file.exists(env_file)) stop("Missing R_environment.txt.")
envs <- read_lines(env_file)

cat("Installing/loading packages from env file:\n")
print(envs)

# Attempt install via pacman (CRAN only)
for (pkg in envs) {
  if (!suppressWarnings(pacman::p_isinstalled(pkg))) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch({
      pacman::p_load(pkg, character.only = TRUE)
    }, error = function(e) {
      cat(sprintf("Trying Bioconductor for %s...\n", pkg))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    })
  } else {
    suppressMessages(library(pkg, character.only = TRUE))
  }
}