# install_r_deps.R

# ── Python: rpy2 ──────────────────────────────────────────────────────────────
message("Installing rpy2...")
system("python -m pip install rpy2 --only-binary=:all:")

# ── R packages ────────────────────────────────────────────────────────────────
lib <- path.expand("~/R/library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", lib=lib, repos="https://cloud.r-project.org/")

BiocManager::install(
    c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb",
      "GenomicRanges", "Biostrings", "DNAshapeR"),
    lib=lib, ask=FALSE, update=FALSE
)

cat("All dependencies installed successfully.\n")
